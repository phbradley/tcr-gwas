library(lme4)
library(data.table)

# First, create a function for regression without observation weighting or varying intercepts

trimming_snp_regression_simple <- function(snps_dataframe, condensed_trimming_dataframe){
    simple_regression_results = data.frame()

    for (snpID in names(snps_dataframe)[-c(1,1185)]){
        sub = NULL
        sub2 = NULL
        regression = NULL
        sub = data.frame(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = as.data.table(merge(sub, condensed_trimming_dataframe, by = "localID"))
        regression = glm(formula = avg_v_gene_trim ~ snp, data = sub2[snp != "NA"])
        simple_regression_results = rbind(simple_regression_results, data.frame(snp = snpID, intercept = as.numeric(coef(regression)[1]), slope = as.numeric(coef(regression)[2]), nonNA_snp_observation_count = as.numeric(nrow(sub2[snp != "NA"]))))
    }

    return(simple_regression_results)
}


# Next, create a function for regression without observation weighting WITH varying intercepts by subject

trimming_snp_regression_by_vgene_subject_varying_intercepts_subject <- function(snps_dataframe, condensed_trimming_dataframe){
    regression_results_subject_varying_int = data.frame()

    for (snpID in names(snps_dataframe)[-c(1,1185)]){
        sub = NULL
        sub2 = NULL
        regression = NULL
        sub = data.frame(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = as.data.table(merge(sub, condensed_trimming_dataframe, by = "localID"))
        regression = lmer(formula = avg_v_gene_trim ~ snp + (1|localID), data = sub2[snp != "NA"])
        regression_results_subject_varying_int = rbind(regression_results_subject_varying_int, data.frame(snp = snpID, intercept = as.numeric(coef(regression)[1]), slope = as.numeric(coef(regression)[2]), nonNA_snp_observation_count = as.numeric(nrow(sub2[snp != "NA"]))))
    }
    
    return(regression_results_subject_varying_int)
}






