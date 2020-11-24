library("lme4")
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

set_regression_formula <- function(trim_type, pca_structure_correction, GENE_CONDITIONING, RANDOM_EFFECTS, by_race = 'False'){
    # define trim type and gene type
    gene_type <<- paste0(substr(trim_type, 1, 1), '_gene')
    
    if (str_split(trim_type, "_")[[1]][2] == 'insert'){
        gene_type1 = paste0(substr(trim_type, 1, 1), '_gene')
        gene_type2 = paste0(substr(trim_type, 2, 2), '_gene')
    }

    # set base regression formula
    if (RANDOM_EFFECTS == 'True'){
        form = formula(get(paste0(trim_type)) ~ snp + (1|localID)) 
    } else {
        form = formula(get(paste0(trim_type)) ~ snp)
    }
    
    if (GENE_CONDITIONING == 'True'){
        if (str_split(trim_type, "_")[[1]][2] == 'insert'){
            form = update(form, ~ . + get(paste0(gene_type1)) + get(paste0(gene_type2)))
        } else {
            form = update(form, ~ . + as.factor(get(paste0(gene_type))))
        }
    } 

    if (pca_structure_correction == 'True'){
       form = update(form, ~ . + EV1 + EV2 + EV3 + EV4 + EV5 + EV6 + EV7 + EV8)
    }

    if (by_race == 'True'){
        form = update(form, ~ . + as.factor(race.g))
    }
    return(form)
}

trimming_regression <- function(snps_dataframe, condensed_trimming_dataframe, CONDENSING, productive, trim_type, REPETITIONS, pca_structure_correction, pca_type, BOOT_CUTOFF, GENE_CONDITIONING, WEIGHTING, RANDOM_EFFECTS, snp_list){
    # set bonferroni correction to us the full group of snps from the gwas (regardless of how many we want to analyze)
    bonferroni = 0.05/35481497

    # subset trimming data to include only productive or not productive entires
    condensed_trimming_dataframe = filter_by_productivity(as.data.frame(condensed_trimming_dataframe), productive)

    # Indicate which snp we are regressing
    snpID = names(snps_dataframe)[-c(1)]
    
    colnames(snps_dataframe) = c('localID', 'snp')
    # merge snp data and trimming data
    snps_trimming_data = merge(snps_dataframe, condensed_trimming_dataframe, by = "localID")
    snps_trimming_data = snps_trimming_data %>% filter(!is.na(snp))

    if (pca_structure_correction != 'False'){
        genotype_pca = read_genotype_pca(pca_type)
        snps_trimming_data = merge(snps_trimming_data, genotype_pca)
    }
    form = set_regression_formula(trim_type, pca_structure_correction, GENE_CONDITIONING, RANDOM_EFFECTS)

    # REGRESSION!
    # remove warning messages (about singularity)
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4), calc.derivs = FALSE)

    if (WEIGHTING == 'True'){
        stopifnot(str_split(trim_type, "_")[[1]][2] == 'trim')
        weight = paste0("weighted_", substr(trim_type, 1, 1), '_gene_count')
        if (RANDOM_EFFECTS == 'True'){
            regression = eval(bquote(lmer(formula = form, data = snps_trimming_data, weights = .(as.name(weight)), control=control)))
        } else {
            regression = eval(bquote(lm(formula = form, data = snps_trimming_data, weights = .(as.name(weight)))))
        }
    } else {
        if (RANDOM_EFFECTS == 'True'){
            regression = lmer(formula = form, data = snps_trimming_data, control=control)
        } else {
            regression = lm(formula = form, data = snps_trimming_data)
        }
    }

    # Calculate slope, intercept 
    if (RANDOM_EFFECTS == 'True'){
        slope = fixef(regression)['snp']
        intercept = fixef(regression)['(Intercept)'] + mean(fixef(regression)[-c(1,2)]) #need to add random effects here...
    } else {
        slope = coef(regression)['snp']
        intercept = coef(regression)['(Intercept)']
    }

    if (slope != 'NA'){
        # Pvalue screen before doing bootstrap (so that we only bootstrap things that may be significant)
        if (CONDENSING == 'by_gene' | CONDENSING == 'by_cdr3') {
            bootstrap_results = bootstrap_screen(regression, RANDOM_EFFECTS)
            if (bootstrap_results$pvalue < BOOT_CUTOFF){
                bootstrap_results = calculate_pvalue(regression, 
                                                     data = snps_trimming_data, 
                                                     cluster_variable = snps_trimming_data$localID, 
                                                     trim_type, 
                                                     RANDOM_EFFECTS, 
                                                     GENE_CONDITIONING, 
                                                     WEIGHTING, 
                                                     REPETITIONS)
            }
        } else {
            bootstrap_results = bootstrap_screen(regression, RANDOM_EFFECTS)
        }
    } else {
        bootstrap_results = data.frame()
    }
    if (substr(snp_list$snp[1], 1, 3) != 'snp'){
        snp_list$snp = paste0('snp', snp_list$snp)
    }
    
    results_temp = merge(snp_list, data.frame(snp = snpID, intercept = intercept, slope = slope), by = 'snp')
    # combine slope, intercept, and pvalue for the specified snp
    regression_results = cbind(results_temp, bootstrap_results)

    return(regression_results)
}
