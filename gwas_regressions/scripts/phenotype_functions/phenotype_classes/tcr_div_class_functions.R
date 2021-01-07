CONDENSING_VARIABLE <<- 'by_subject'
CONDITIONING_VARIABLE <<- FALSE
INFER_MISSING_D_GENE <<- 'False'
REPETITIONS <<- FALSE
BOOTSTRAP_PVALUE_CUTOFF <<- FALSE
PCA_COUNT <<- 8
 

condense_individual_tcr_repertoire_data <- function(localID, tcr_repertoire_dataframe){
    tcr_repertoire = data.table(localID = localID)
    tcr_repertoire$tcr_count = nrow(tcr_repertoire_dataframe)

    tcr_div_file = fread(TCR_DIV_FILE)
    setnames(tcr_div_file, 'hipid', 'localID')
    setnames(tcr_div_file, 'tcrdiv', 'tcr_div')
    tcr_repertoire = merge(tcr_repertoire, tcr_div_file, by = 'localID')

    return(tcr_repertoire)
}

condense_all_tcr_repertoire_data <- function(){
    files = list.files(TCR_REPERTOIRE_DATA_DIRECTORY, pattern = "*.tsv", full.names=TRUE)
    tcr_rep_data = data.table()        

    count = 0

    registerDoParallel(cores=NCPU)
    tcr_rep_data = foreach(file = files, .combine = 'rbind') %dopar% {
        count = count + 1
        file_data = fread(file)
        localID = extract_subject_ID(file)
        print(paste0('processing ', count, ' of ', length(files)))
        condense_individual_tcr_repertoire_data(localID, file_data)
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
    if (length(unique(regression_data[[snp]][!is.na(regression_data[[snp]])])) <= 1){
        next
    }
    regression = eval(bquote(lm(formula = formula, data = regression_data)))
    regression_results = calculate_regression_results(snp, regression, regression_data)
    regression_results_together = rbind(regression_results_together, regression_results)
    print(paste0('finished regression for snp ', index, ' of ', length(snps), ' for ', PHENOTYPE))
    return(regression_results_together)    
}


