library("lme4")
library("data.table")

source("trimming_bootstrap_functions.R")

# Weighted mixed model function

trimming_snp_regression_weighted_varying_int_subject <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, repetitions, bonferroni){
    # remove warning messages (about singularity)
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
    # set bonferroni correction to us the full group of snps from the gwas (regardless of how many we want to analyze)
    bonferroni = 0.05/35481497

    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    simple_regression_results = data.table()
    bootstrap_results = data.table()

    # subset trimming data to include only productive or not productive entires
    if (trim_type == "v_trim" | trim_type == "j_trim" | trim_type == "vj_insert"){
        if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
        } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
        }
    } else if (trim_type == "d0_trim" | trim_type == "d1_trim" | trim_type == "dj_insert" | trim_type == "vd_insert"){
        if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "TRUE"]
        } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "FALSE"]
        }
    } 
    
 
    # For each snpID:
    for (snpID in names(snps_dataframe)[-c(1)]){
        sub = NULL
        sub2 = NULL
        regression = NULL

        # merge snp data and trimming data
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    
        # set regression formula given 
        form = formula(get(paste0(trim_type)) ~ snp + (1|localID))

        # set weights for regression
        if (trim_type =='d1_trim' | trim_type =='d0_trim'){
            weight = paste0("weighted_d_gene_count")
        } else {
            weight = paste0("weighted_", str_split(trim_type, "_")[[1]][1], "_gene_count")
        }

        # REGRESSION!
        regression = lmer(formula = form, data = sub2[snp != "NA"], weights = get(weight), control=control)
        
        # Calculate slope, intercept 
        intercept = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][1]))))
        slope = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][2]))))

        if (slope != "NA"){
            bootstrap_results = calculate_pvalue(regression, data = sub2[snp != "NA"], cluster_variable = sub2[snp != "NA"]$localID, repetitions)
            if (bootstrap_results[2]<bonferroni){
                bootstrap_results = calculate_pvalue(regression, data = sub2[snp != "NA"], cluster_variable = sub2[snp != "NA"]$localID, repetitions=10)
            } else {
                bootstrap_results = bootstrap_results
            }
        } else {
            bootstrap_results = data.table()
        }
        together = cbind(data.table(snp = snpID, intercept = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][1])))), slope = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][2]))))), bootstrap_results)
        
        # Combine snpID, intercept, slope, etc.
        results = rbind(simple_regression_results, together)
    }
    return(results)
}