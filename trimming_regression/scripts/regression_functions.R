library("lme4")
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

# this script does a lmer regression including fixed and random effects to condition out the effects mediated by gene choice

trimming_regression <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, gene_type, bootstrap_repetitions, pca_structure_correction, gene_conditioning, weighting, random_effects, snp_list){
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

    genotype_pca = read_genotype_pca()
    snps_trimming_data = merge(snps_trimming_data, genotype_pca)

    # set regression formula
    # set base formula
    if (random_effects == 'True'){
        form = formula(get(paste0(trim_type)) ~ snp + (1|localID))
    } else {
        form = formula(get(paste0(trim_type)) ~ snp)
    }
    
    if (gene_conditioning == 'True'){
        if (str_split(trim_type, "_")[[1]][2] == 'insert'){
            form = update(form, ~ . + get(paste0(gene_type1)) + get(paste0(gene_type2)))
        } else {
            form = update(form, ~ . + get(paste0(gene_type)))
        }
    } else {
        weight = 'tcr_count'
    }

    if (pca_structure_correction == 'True'){
       form = update(form, ~ . + EV1 + EV2 + EV3 + EV4 + EV5 + EV6 + EV7 + EV8 + EV9 + EV10)
    } 

    # REGRESSION!
    # remove warning messages (about singularity)
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4), calc.derivs = FALSE)

    if (weighting == 'True'){
        if (random_effects == 'True'){
            regression = lmer(formula = form, data = snps_trimming_data, weights = get(weight), control=control)
        } else {
            regression = lm(formula = form, data = snps_trimming_data, weights = get(weight))
        }
    } else {
        if (random_effects == 'True'){
            regression = lmer(formula = form, data = snps_trimming_data, control=control)
        } else {
            regression = lm(formula = form, data = snps_trimming_data)
        }
    }

    # Calculate slope, intercept 
    if (random_effects == 'True'){
        slope = fixef(regression)['snp']
        intercept = fixef(regression)['(Intercept)'] + mean(fixef(regression)[-c(1,2)]) #need to add random effects here...
    } else {
        slope = coef(regression)['snp']
        intercept = coef(regression)['(Intercept)']
    }

    if (slope != 'NA'){
        # Pvalue screen before doing bootstrap (so that we only bootstrap things that may be significant)
        if (bootstrap_repetitions != 0){
            bootstrap_results = calculate_pvalue(regression, data = snps_trimming_data, cluster_variable = snps_trimming_data$localID, trim_type, varying_int = random_effects, gene_conditioning, weighting, repetitions = bootstrap_repetitions)
        } else {
            bootstrap_results = bootstrap_screen(regression, random_effects)
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
