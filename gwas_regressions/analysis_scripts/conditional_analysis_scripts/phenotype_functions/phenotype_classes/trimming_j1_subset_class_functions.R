ALLELE_STATUS_CORRECTION <<- NA

condense_individual_tcr_repertoire_data <- function(tcr_repertoire_dataframe){
    tcr_repertoire_dataframe[, productivity_tcr_count := .N, by = .(localID, productive)]
    tcr_repertoire_dataframe = tcr_repertoire_dataframe[substring(j_gene, 1, 5) == 'TRBJ1']
    
    tcr_repertoire_dataframe$tcr_count = nrow(tcr_repertoire_dataframe)

    names = c('localID', 'productive')
    
    condensed_tcr_repertoire_data = tcr_repertoire_dataframe[, lapply(.SD, mean), by = mget(names), .SDcols = sapply(tcr_repertoire_dataframe, is.numeric)]

    setnames(condensed_tcr_repertoire_data, gsub('_j1_subset', '', PHENOTYPE), PHENOTYPE)
    valid_columns = c(names, PHENOTYPE, 'tcr_count', 'productivity_tcr_count')

    return(condensed_tcr_repertoire_data[,..valid_columns])
}

condense_all_tcr_repertoire_data <- function(){
    files = list.files(TCR_REPERTOIRE_DATA_DIRECTORY, pattern = "*.tsv", full.names=TRUE)
    tcr_rep_data = data.table()        

    count = 0

    registerDoParallel(cores=NCPU)
    tcr_rep_data = foreach(file = files, .combine = 'rbind') %dopar% {
        count = count + 1
        file_data = fread(file)
        file_data$localID = extract_subject_ID(file)
        if (KEEP_MISSING_D_GENE == 'False'){
            file_data = file_data[d_gene != '-']
        }
        print(paste0('processing ', count, ' of ', length(files)))
        condense_individual_tcr_repertoire_data(file_data)
    }
    stopImplicitCluster()
    filename = generate_condensed_tcr_repertoire_file_name()
    write.table(tcr_rep_data, file = filename, quote=FALSE, sep='\t', col.names = NA)
}

set_regression_formula <- function(snp, conditional_snp_list = NULL){
    pca_covariates = paste0('EV', seq(1, PCA_COUNT), collapse = '+')
    formula_string = paste0(PHENOTYPE, ' ~ `', snp, '`+')
    if (!is.null(conditional_snp_list)){
        formula_string = paste0(formula_string, '`', paste0(conditional_snp_list, collapse = '`+`'), '`+')
    }
    formula = formula(paste(formula_string, pca_covariates, sep = ''))
    return(formula)
}


calculate_regression_results <- function(snp, regression, regression_data_filtered_by_productivity, conditional_snp_list){
    regression_results = data.table(snp = snp, 
                                    slope = summary(regression)$coefficients[paste0('`', snp, '`'),'Estimate'], 
                                    standard_error = summary(regression)$coefficients[paste0('`', snp, '`'),'Std. Error'], 
                                    pvalue = summary(regression)$coefficients[paste0('`', snp, '`'),'Pr(>|t|)'], 
                                    parameter = 'snp', 
                                    phenotype = PHENOTYPE, 
                                    bootstraps = 0)
    for (cond_snp in seq(1, length(conditional_snp_list))){
        regression_results[, paste0('conditional_snp_', cond_snp) := conditional_snp_list[cond_snp]]
        regression_results[, paste0('conditional_snp_', cond_snp, '_slope') := summary(regression)$coefficients[paste0('`', conditional_snp_list[cond_snp], '`'),'Estimate']]
        regression_results[, paste0('conditional_snp_', cond_snp, '_standard_error') := summary(regression)$coefficients[paste0('`', conditional_snp_list[cond_snp], '`'),'Std. Error']]
        regression_results[, paste0('conditional_snp_', cond_snp, '_pvalue') := summary(regression)$coefficients[paste0('`', conditional_snp_list[cond_snp], '`'),'Pr(>|t|)']]
    }
    return(regression_results)
}

   
conditional_regress <- function(snp, snps, regression_data, conditional_snp_list, productivity){
    index = which(snps == snp)
    formula = set_regression_formula(snp, conditional_snp_list)
    regression_results_together = data.table()
    all_snps = c(snp, conditional_snp_list)
    non_na_data = regression_data[productive == productivity][, ..all_snps][complete.cases(regression_data[productive == productivity][, ..all_snps]),]
    if (nrow(non_na_data) > 0 & length(unique(as.list(non_na_data))) == length(all_snps)){
        regression = eval(bquote(lm(formula = formula, data = regression_data[productive == productivity])))
        regression_results = calculate_regression_results(snp, regression, regression_data[productive == productivity], conditional_snp_list)
        regression_results$productive = productivity
        regression_results$conditional_snps = paste0(conditional_snp_list, collapse = ', ')
        regression_results_together = rbind(regression_results_together, regression_results)
        print(paste0('finished regression for snp ', index, ' of ', length(snps), ' for ', PHENOTYPE))
    } else {
        print(paste0('NO regression possible for snp ', index, ' of ', length(snps), ' for ', PHENOTYPE)) 
    }
    return(regression_results_together)    
}
