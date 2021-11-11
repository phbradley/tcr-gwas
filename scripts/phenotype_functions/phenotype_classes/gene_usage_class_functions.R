ALLELE_STATUS_CORRECTION <<- NA

get_gene_from_gene_allele <- function(gene_allele_list){
    split_gene_allele = strsplit(gene_allele_list, '*', fixed = TRUE)
    genes = sapply(split_gene_allele, function(x) x[1])
    return(genes)
}

condense_individual_tcr_repertoire_data <- function(tcr_repertoire_dataframe){
    tcr_repertoire_dataframe$tcr_count = nrow(tcr_repertoire_dataframe)
    tcr_repertoire_dataframe[, productivity_tcr_count := .N, by = .(localID, productive)]

    gene_usage = data.table()

    for (gene_type in c('v_gene', 'd_gene', 'j_gene')){
        if (gene_type == 'd_gene' & KEEP_MISSING_D_GENE == 'True'){
            temp_tcr_repertoire_dataframe = infer_d_gene(tcr_repertoire_dataframe)
        } else if (gene_type == 'd_gene' & KEEP_MISSING_D_GENE == 'False'){
            temp_tcr_repertoire_dataframe = tcr_repertoire_dataframe[d_gene != '-']
        } else {
            temp_tcr_repertoire_dataframe = tcr_repertoire_dataframe
        }

        temp_tcr_repertoire_dataframe[,gene:=get_gene_from_gene_allele(get(gene_type))]
        gene_freqs = temp_tcr_repertoire_dataframe[,.N, by = .(gene, localID, productive, tcr_count, productivity_tcr_count)]
        colnames(gene_freqs) = c('gene', 'localID', 'productive', 'tcr_count', 'productivity_tcr_count', 'gene_count')
        gene_freqs[, gene_usage := gene_count/tcr_count]
        gene_usage = rbind(gene_usage, gene_freqs)
    }
    # remove orphan genes
    gene_usage = gene_usage[!(gene %like% 'OR')]
    return(gene_usage)
}

construct_complete_dataframe <- function(subjects, genes){
    together = data.table()
    genes_dt = data.table(gene = genes)
    for (subject in subjects){
        for (prod in c(TRUE, FALSE)){
            genes_dt$productive = prod
            genes_dt$localID = subject
            together = rbind(together, genes_dt)
        }
    }
    return(together)
}

fill_in_zero_frequencies <- function(compiled_dataframe){
    meta_subject_data = unique(compiled_dataframe[, c('localID', 'tcr_count', 'productivity_tcr_count', 'productive')])
    all_genes = unique(compiled_dataframe$gene)
    all_subjects = unique(compiled_dataframe$localID)
    complete_dataframe = construct_complete_dataframe(all_subjects, all_genes)
    complete_dataframe = merge(complete_dataframe, meta_subject_data, by = c('localID', 'productive'))
    complete_dataframe = merge(complete_dataframe, compiled_dataframe, by = colnames(complete_dataframe), all.x = TRUE)
    complete_dataframe[is.na(gene_count), gene_count := 0]
    complete_dataframe[is.na(gene_usage), gene_usage := 0]
    return(complete_dataframe)
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
    tcr_rep_data = fill_in_zero_frequencies(tcr_rep_data)
    filename = generate_condensed_tcr_repertoire_file_name()
    write.table(tcr_rep_data, file = filename, quote=FALSE, sep='\t', col.names = NA)
}

set_regression_formula <- function(snp){
    if (!is.na(PCA_COUNT)){
        pca_covariates = paste0('EV', seq(1, PCA_COUNT), collapse = '+')
        formula_string = paste0(PHENOTYPE, ' ~ `', snp, '`+')
        formula = formula(paste(formula_string, pca_covariates, sep = ''))
    } else {
        formula = formula(paste0(PHENOTYPE, ' ~ `', snp, '`'))
    }
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
    if (PHENOTYPE != 'gene_usage'){
        setnames(regression_data, 'gene_usage', PHENOTYPE, skip_absent = TRUE)
    }
    for (productivity in c('TRUE', 'FALSE')){
        for(gene_name in unique(regression_data$gene)){
            # if all subjects have the same genotype, skip
            if (length(unique(regression_data[productive == productivity & gene == gene_name][[snp]][!is.na(regression_data[productive == productivity & gene == gene_name][[snp]])])) <= 1){
                next
            }
            # if the number of data points is less than or equal to the number of parameters to be fit, skip
            if (length(regression_data[productive == productivity & gene == gene_name][[snp]][!is.na(regression_data[productive == productivity & gene == gene_name][[snp]])]) <= 10){
                next
            }
            regression = eval(bquote(lm(formula = formula, data = regression_data[productive == productivity & gene == gene_name])))
            regression_results = calculate_regression_results(snp, regression, regression_data[productive == productivity & gene == gene_name])
            regression_results$productive = productivity
            regression_results$gene = gene_name
            regression_results_together = rbind(regression_results_together, regression_results)
        }
    }
    print(paste0('finished regression for snp ', index, ' of ', length(snps), ' for ', PHENOTYPE))
    return(regression_results_together)    
}


