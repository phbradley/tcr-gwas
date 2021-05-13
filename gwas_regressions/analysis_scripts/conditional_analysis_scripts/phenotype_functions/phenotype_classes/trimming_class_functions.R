ALLELE_STATUS_CORRECTION <<- NA

condense_individual_tcr_repertoire_data <- function(tcr_repertoire_dataframe){
    cdr3 = combine_genes_by_common_cdr3()
    tcr_repertoire_data = merge(tcr_repertoire_dataframe, cdr3, by.x = GENE_TYPE, by.y = 'id')
    tcr_repertoire_data$tcr_count = nrow(tcr_repertoire_data)
    tcr_repertoire_data[, productivity_tcr_count := .N, by = .(localID, productive)]

    names = c('localID', paste(CONDITIONING_VARIABLE), 'productive')
    
    condensed_tcr_repertoire_data = tcr_repertoire_data[, paste0(GENE_TYPE, '_count') := .N, by = mget(names)][, lapply(.SD, mean), by = mget(names), .SDcols = sapply(tcr_repertoire_data, is.numeric)]

    condensed_tcr_repertoire_data[[paste0('weighted_', GENE_TYPE, '_count')]] = condensed_tcr_repertoire_data[[paste0(GENE_TYPE, '_count')]]/nrow(tcr_repertoire_data)

    valid_columns = c(names, PHENOTYPE, 'tcr_count', 'productivity_tcr_count', paste0(GENE_TYPE, '_count'), paste0('weighted_', GENE_TYPE, '_count'))

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
        if (KEEP_MISSING_D_GENE == 'False' & CHAIN == 'beta'){
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
    formula_string = paste0(PHENOTYPE, ' ~ `', snp, '` + ', 'as.factor(', CONDITIONING_VARIABLE, ') +')
    if (!is.null(conditional_snp_list)){
        formula_string = paste0(formula_string, '`', paste0(conditional_snp_list, collapse = '`+`'), '`+')
    }
    formula = formula(paste(formula_string, pca_covariates, sep = ''))
    return(formula)
}

bootstrap_by_localID <- function(snp, list_of_snps_in_regression, regression, regression_data_filtered_by_productivity){
    subjects = names(table(regression_data_filtered_by_productivity$localID))
    bootstrap_results = matrix(NA, nrow=REPETITIONS, ncol = 1)
    weight = paste0('weighted_', GENE_TYPE, '_count')

    for (i in 1:REPETITIONS){
        sample_subjects = sample(1:length(subjects), length(subjects), replace = TRUE)
        subjects_by_sampling = subjects[sample_subjects]
        contingency_table_samples = table(subjects_by_sampling)
        bootstrap_data = NULL
        for (subject_sample_count in 1:max(contingency_table_samples)){
            data_subset = regression_data_filtered_by_productivity[regression_data_filtered_by_productivity$localID %in% names(contingency_table_samples[contingency_table_samples %in% subject_sample_count])]
            for (count in 1:subject_sample_count){
                bootstrap_data = rbind(bootstrap_data, data_subset)
            }
        }
        
        if (length(unique(bootstrap_data[[snp]])[!is.na(unique(bootstrap_data[[snp]]))])<= 1){
            next
        }
        
        all_snps = list_of_snps_in_regression     
        non_na_data = bootstrap_data[, ..all_snps][complete.cases(bootstrap_data[, ..all_snps]),]
        if (nrow(non_na_data) == 0 | length(unique(as.list(non_na_data))) == 1){
            next
        }

        unique_genotypes_per_snp = lengths(lapply(non_na_data, unique))
        nonsingular_unique_genotypes_per_snp = unique_genotypes_per_snp[unique_genotypes_per_snp != 1]
        if (length(nonsingular_unique_genotypes_per_snp) != length(unique_genotypes_per_snp)){
            next
        }

        formula = formula(regression)

        bootstrap_regression = eval(bquote(lm(formula = formula, data = bootstrap_data, weights = .(as.name(weight)))))
        bootstrap_results[i,] = summary(bootstrap_regression)$coefficients[paste0('`', snp, '`'),'Estimate']
    }
    standard_error = apply(na.omit(bootstrap_results),2,sd)
    return(standard_error)
}

recalculate_pvalue <- function(slope, standard_error){
    zscore = slope/standard_error
    pvalue = 2*pnorm(-abs(zscore))
    return(pvalue)
}

calculate_regression_results <- function(snp, regression, regression_data_filtered_by_productivity, conditional_snp_list){
    regression_results = data.table(snp = snp, 
                                    slope = summary(regression)$coefficients[paste0('`', snp, '`'),'Estimate'], 
                                    standard_error = summary(regression)$coefficients[paste0('`', snp, '`'),'Std. Error'], 
                                    pvalue = summary(regression)$coefficients[paste0('`', snp, '`'),'Pr(>|t|)'], 
                                    parameter = 'snp', 
                                    phenotype = PHENOTYPE, 
                                    bootstraps = 0)
    if (regression_results$pvalue < BOOTSTRAP_PVALUE_CUTOFF){
        regression_results$standard_error = bootstrap_by_localID(snp, c(snp, conditional_snp_list), regression, regression_data_filtered_by_productivity)        
        regression_results$pvalue = recalculate_pvalue(regression_results$slope, regression_results$standard_error)
        regression_results$bootstraps = REPETITIONS
    }
    for (cond_snp in seq(1, length(conditional_snp_list))){
        regression_results[, paste0('conditional_snp_', cond_snp) := conditional_snp_list[cond_snp]]
        regression_results[, paste0('conditional_snp_', cond_snp, '_slope') := summary(regression)$coefficients[paste0('`', conditional_snp_list[cond_snp], '`'),'Estimate']]
        regression_results[, paste0('conditional_snp_', cond_snp, '_standard_error') := summary(regression)$coefficients[paste0('`', conditional_snp_list[cond_snp], '`'),'Std. Error']]
        regression_results[, paste0('conditional_snp_', cond_snp, '_pvalue') := summary(regression)$coefficients[paste0('`', conditional_snp_list[cond_snp], '`'),'Pr(>|t|)']]

        if (regression_results[[paste0('conditional_snp_', cond_snp, '_pvalue')]] < BOOTSTRAP_PVALUE_CUTOFF){
            regression_results[[paste0('conditional_snp_', cond_snp, '_standard_error')]] = bootstrap_by_localID(paste0(conditional_snp_list[cond_snp]), c(snp, conditional_snp_list), regression, regression_data_filtered_by_productivity)
            regression_results[[paste0('conditional_snp_', cond_snp, '_pvalue')]] = recalculate_pvalue(regression_results[[paste0('conditional_snp_', cond_snp, '_slope')]], regression_results[[paste0('conditional_snp_', cond_snp, '_standard_error')]])
        }
    }
    return(regression_results)
}

conditional_regress <- function(snp, snps, regression_data, conditional_snp_list, productivity){
    index = which(snps == snp)
    formula = set_regression_formula(snp, conditional_snp_list)
    weight = paste0('weighted_', GENE_TYPE, '_count')
    regression_results_together = data.table()
    all_snps = c(snp, conditional_snp_list)
    non_na_data = regression_data[productive == productivity][, ..all_snps][complete.cases(regression_data[productive == productivity][, ..all_snps]),]
    if (nrow(non_na_data) > 0 & length(unique(as.list(non_na_data))) > 1){
        regression = eval(bquote(lm(formula = formula, data = regression_data[productive == productivity], weights = .(as.name(weight)))))
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

