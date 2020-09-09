library("lme4")

source("/home/mrussel2/tcr-gwas/trimming_regression/bootstrap_functions.R")

trimming_regression <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, gene_type, bootstrap_repetitions, gene_conditioning, weighting, snp_list){
    # set bonferroni correction to us the full group of snps from the gwas (regardless of how many we want to analyze)
    bonferroni = 0.05/35481497

    # subset trimming data to include only productive or not productive entires
    condensed_trimming_dataframe = filter_by_productivity(as.data.frame(condensed_trimming_dataframe), productive)

    # Indicate which snp we are regressing
    snpID = names(snps_dataframe)[-c(1)]

    # define trim type and gene type
    if (gene_type == 'same'){
        gene_type = paste0(substr(trim_type, 1, 1), '_gene')
        weight = paste0("weighted_", substr(trim_type, 1, 1), '_gene_count')
    }
    gene_type = paste0(gene_type)
    if (str_split(trim_type, "_")[[1]][2] == 'insert'){
        gene_type1 = paste0(substr(trim_type, 1, 1), '_gene')
        gene_type2 = paste0(substr(trim_type, 2, 2), '_gene')
        weight = paste0("weighted_", substr(trim_type, 1, 2), '_gene_count')
    }

    colnames(snps_dataframe) = c('localID', 'snp')
    # merge snp data and trimming data
    snps_trimming_data = merge(snps_dataframe, condensed_trimming_dataframe, by = "localID")
    snps_trimming_data = snps_trimming_data %>% filter(!is.na(snp))

    # set regression formula
    # set base formula
    form = formula(get(paste0(trim_type)) ~ snp + (1|localID))
    
    
    if (gene_conditioning == 'True'){
        if (str_split(trim_type, "_")[[1]][2] == 'insert'){
            form = update(form, ~ . + get(paste0(gene_type1)) + get(paste0(gene_type2)))
        } else {
            form = update(form, ~ . + get(paste0(gene_type)))
        }
    } 

    # REGRESSION!
    # remove warning messages (about singularity)
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4), calc.derivs = FALSE)

    if (weighting == 'True'){
        regression = lmer(formula = form, data = snps_trimming_data, weights = get(weight), control=control)
    } else {
        regression = lmer(formula = form, data = snps_trimming_data, control=control)
    }

    # Calculate slope, intercept 
    # Add the Intercept term with a mean of the gene specific intercept
    intercept = summary(regression)$coefficients[,'Estimate']['(Intercept)'] + mean(summary(regression)$coefficients[,'Estimate'][-c(1,2)])
    slope = summary(regression)$coefficients[,'Estimate']['snp']

    if (slope != 'NA'){
        # Pvalue screen before doing bootstrap (so that we only bootstrap things that may be significant)
        if (bootstrap_repetitions != 0){
            bootstrap_results = calculate_pvalue(regression, data = snps_trimming_data[snp != "NA"], cluster_variable = snps_trimming_data[snp != "NA"]$localID, trim_type, varying_int = 'True', weighting, repetitions = bootstrap_repetitions)
        } else {
            bootstrap_results = bootstrap_screen(regression)
        }
        #boot_screen = bootstrap_screen(regression)
        #if (bootstrap_repetitions != 0 & boot_screen[2]< (bonferroni)){
        #    bootstrap_results = calculate_pvalue(regression, data = snps_trimming_data[snp != "NA"], cluster_variable = snps_trimming_data[snp != "NA"]$localID, trim_type, varying_int, weighting, bootstrap_repetitions)
        #} else {
        #    bootstrap_results = boot_screen
        #}
    } else {
        bootstrap_results = data.frame()
    }

    snp_list$snp = paste0('snp', snp_list$snp)
    results_temp = merge(snp_list, data.frame(snp = snpID, intercept = intercept, slope = slope), by = 'snp')
    # combine slope, intercept, and pvalue for the specified snp
    regression_results = cbind(results_temp, bootstrap_results)

    return(regression_results)
}


generate_file_name <- function(snp_id_list, trim_type, gene_type, productivity, gene_conditioning, weighting, condensing, repetitions){
    prod = ifelse(productivity == 'True', 'productive', 'NOT_productive')
    gene = ifelse(gene_conditioning == 'True', 'with_gene', '')
    weight = ifelse(weighting == 'True', '_with_weighting', '')
    if ((substr(gene_type, 1, 1) == substr(trim_type, 1, 1)) | (substr(trim_type, 4, 5) == 'in')){
        name = paste0('regression_bootstrap_results/', prod, '/', trim_type, '/', trim_type, '_', prod, '_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_lmer_', gene, weight, '_condensing_', condensing, '_', repetitions, '_bootstraps.tsv') 
    } else {
        name = paste0('regression_bootstrap_results/', prod, '/crosses/', trim_type, '_', gene_type, '_', prod, '_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_lmer_', gene, weight, '_condensing_', condensing, '_', repetitions, '_bootstraps.tsv') 
    }
    return(name)
}
        
    
    
generate_file_name_no_reps <- function(snp_id_list, trim_type, productivity, gene_conditioning, weighting, condensing){
    prod = ifelse(productivity == 'True', 'productive', 'NOT_productive')
    gene = ifelse(gene_conditioning == 'True', 'with_gene', '')
    weight = ifelse(weighting == 'True', '_with_weighting', '')

    name = paste0('regression_bootstrap_results/', prod, '/', trim_type, '/', trim_type, '_', prod, '_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_lmer_', gene, weight, '_condensing_', condensing, '.tsv') 

    return(name)
}