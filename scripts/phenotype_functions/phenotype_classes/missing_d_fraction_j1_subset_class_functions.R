ALLELE_STATUS_CORRECTION <<- NA

condense_individual_tcr_repertoire_data <- function(tcr_repertoire_dataframe){
    tcr_repertoire_dataframe[, productivity_tcr_count := .N, by = .(localID, productive)]
    tcr_repertoire_dataframe = tcr_repertoire_dataframe[substring(j_gene, 1, 5) == 'TRBJ1']
    
    tcr_repertoire_dataframe$tcr_count = nrow(tcr_repertoire_dataframe)
    tcr_repertoire_dataframe[, productivity_tcr_count_j1_subset := .N, by = .(localID, productive)]

    names = c('localID', 'productive')
    
    condensed_tcr_repertoire_data = tcr_repertoire_dataframe[d_gene == '-', .N, by = .(localID, productive, productivity_tcr_count, tcr_count,  productivity_tcr_count_j1_subset)]
    
    setnames(condensed_tcr_repertoire_data, 'N', 'missing_d_count_j1_subset')
    condensed_tcr_repertoire_data[, missing_d_fraction_j1_subset := missing_d_count_j1_subset/productivity_tcr_count_j1_subset]

    return(condensed_tcr_repertoire_data)
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
        print(paste0('processing ', count, ' of ', length(files)))
        condense_individual_tcr_repertoire_data(file_data)
    }
    stopImplicitCluster()
    filename = generate_condensed_tcr_repertoire_file_name()
    write.table(tcr_rep_data, file = filename, quote=FALSE, sep='\t', col.names = NA)
}

set_regression_formula <- function(snp){
    pca_covariates = paste0('EV', seq(1, PCA_COUNT), collapse = '+')
    formula_string = paste0(PHENOTYPE, ' ~ `', snp, '`+')
    formula = formula(paste(formula_string, pca_covariates, sep = ''))
    return(formula)
}


calculate_regression_results <- function(snp, regression, regression_data_filtered_by_productivity){
    regression_results = data.table(snp = snp, 
                                    slope = summary(regression)$coefficients[paste0('`', snp, '`'),'Estimate'], 
                                    standard_error = summary(regression)$coefficients[paste0('`', snp, '`'),'Std. Error'], 
                                    pvalue = summary(regression)$coefficients[paste0('`', snp, '`'),'Pr(>|t|)'], 
                                    parameter = 'snp', 
                                    phenotype = PHENOTYPE, 
                                    bootstraps = 0)
    return(regression_results)
}
    
regress <- function(snp, snps, regression_data){
    index = which(snps == snp)
    formula = set_regression_formula(snp)
    regression_results_together = data.table()
    for (productivity in c('TRUE', 'FALSE')){
        if (length(unique(regression_data[productive == productivity][[snp]][!is.na(regression_data[productive == productivity][[snp]])])) <= 1){
            next
        }
        regression = eval(bquote(lm(formula = formula, data = regression_data[productive == productivity])))
        regression_results = calculate_regression_results(snp, regression, regression_data[productive == productivity])
        regression_results$productive = productivity
        regression_results_together = rbind(regression_results_together, regression_results)
    }
    print(paste0('finished regression for snp ', index, ' of ', length(snps), ' for ', PHENOTYPE))
    return(regression_results_together)    
}



