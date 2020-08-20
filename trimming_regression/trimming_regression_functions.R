library("lme4")
library("data.table")

source("trimming_bootstrap_functions.R")

trimming_regression <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, bootstrap_repetitions, gene_conditioning, weighting){
    # Indicate that we want to allow for varying intercepts (i.e. random effects by individual)
    varying_int = "True"
    # set bonferroni correction to us the full group of snps from the gwas (regardless of how many we want to analyze)
    bonferroni = 0.05/35481497

    # subset trimming data to include only productive or not productive entires
    if (productive == "True"){
        condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)[productive == "TRUE"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)[productive == "FALSE"]
    }

    # Indicate which snp we are regressing
    snpID = names(snps_dataframe)[-c(1)]

    if (trim_type =='d1_trim' | trim_type =='d0_trim'){
        weight = paste0("weighted_d_gene_count")
        gene_type = paste0("d_gene")
    } else {
        weight = paste0("weighted_", str_split(trim_type, "_")[[1]][1], "_gene_count")
        gene_type = paste0(str_split(trim_type, "_")[[1]][1], "_gene")
        if (str_split(trim_type, "_")[[1]][2] == 'insert'){
            gene_type1 = paste0(substr(trim_type, 1, 1), '_gene')
            gene_type2 = paste0(substr(trim_type, 2, 2), '_gene')
        }
    }

    colnames(snps_dataframe) = c('localID', 'snp')
    # merge snp data and trimming data
    snps_trimming_data = as.data.table(merge(snps_dataframe, condensed_trimming_dataframe, by = "localID"))
    snps_trimming_data = snps_trimming_data[snp != 'NA']

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
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))

    if (weighting == 'True'){
        regression = lmer(formula = form, data = snps_trimming_data[snp != "NA"], weights = get(weight), control=control)
    } else {
        regression = lmer(formula = form, data = snps_trimming_data[snp != "NA"], control=control)
    }

    # Calculate slope, intercept 
    # Add the Intercept term with a mean of the gene specific intercept
    intercept = summary(regression)$coefficients[,'Estimate']['(Intercept)'] + mean(summary(regression)$coefficients[,'Estimate'][-c(1,2)])
    slope = summary(regression)$coefficients[,'Estimate']['snp']

    if (slope != 'NA'){
        # Pvalue screen before doing bootstrap (so that we only bootstrap things that may be significant)
        boot_screen = bootstrap_screen(regression)
        if (bootstrap_repetitions != 0 & boot_screen[2]< (bonferroni)){
            bootstrap_results = calculate_pvalue(regression, data = snps_trimming_data[snp != "NA"], cluster_variable = snps_trimming_data[snp != "NA"]$localID, trim_type, varying_int, weighting, bootstrap_repetitions)
        } else {
            bootstrap_results = boot_screen
        }
    } else {
        bootstrap_results = data.table()
    }

    # combine slope, intercept, and pvalue for the specified snp
    regression_results = cbind(data.table(snp = snpID, intercept = intercept, slope = slope), bootstrap_results)

    return(regression_results)
}


generate_file_name <- function(snp_id_list, trim_type, productivity, gene_conditioning, weighting, condensing, repetitions){
    prod = ifelse(productivity == 'True', 'productive', 'NOT_productive')
    gene = ifelse(gene_conditioning == 'True', 'with_gene', '')
    weight = ifelse(weighting == 'True', '_with_weighting', '')

    name = paste0('regression_bootstrap_results/', prod, '/', trim_type, '/', trim_type, '_', prod, '_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_lmer_', gene, weight, '_condensing_', condensing, '_', repetitions, '_bootstraps.tsv') 

    return(name)
}
        
    
    
