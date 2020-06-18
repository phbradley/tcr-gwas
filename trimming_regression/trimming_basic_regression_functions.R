library(lme4)
library(data.table)

source("trimming_bootstrap_functions.R")
# First, create a function for regression without observation weighting or varying intercepts

trimming_snp_regression_simple <- function(snps_dataframe, condensed_trimming_dataframe, productive){
    simple_regression_results = data.frame()
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    for (snpID in names(snps_dataframe)[-c(1,1185)]){
        sub = NULL
        sub2 = NULL
        regression = NULL
        sub = data.frame(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = as.data.table(merge(sub, condensed_trimming_dataframe, by = "localID"))
        regression = glm(formula = v_trim ~ snp, data = sub2[snp != "NA"])
        simple_regression_results = rbind(simple_regression_results, data.frame(snp = snpID, intercept = as.numeric(coef(regression)[1]), slope = as.numeric(coef(regression)[2]), nonNA_snp_observation_count = as.numeric(nrow(sub2[snp != "NA"]))))
    }

    return(simple_regression_results)
}


# Next, create a function for regression without observation weighting WITH varying intercepts by subject

trimming_snp_regression_by_vgene_subject_varying_intercepts_subject <- function(snps_dataframe, condensed_trimming_dataframe, productive){
    regression_results_subject_varying_int = data.frame()
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    for (snpID in names(snps_dataframe)[-c(1,1185)]){
        sub = NULL
        sub2 = NULL
        regression = NULL
        sub = data.frame(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = as.data.table(merge(sub, condensed_trimming_dataframe, by = "localID"))
        regression = lmer(formula = v_trim ~ snp + (1|localID), data = sub2[snp != "NA"])
        regression_results_subject_varying_int = rbind(regression_results_subject_varying_int, data.frame(snp = snpID, intercept = as.numeric(coef(regression)[1]), slope = as.numeric(coef(regression)[2]), nonNA_snp_observation_count = as.numeric(nrow(sub2[snp != "NA"]))))
    }
    
    return(regression_results_subject_varying_int)
}


trimming_snp_regression_weighted <- function(snps_dataframe, condensed_trimming_dataframe, productive, gene_type, repetitions){
    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    simple_regression_results = data.table()
    bootstrap_results = data.table()
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    for (snpID in names(snps_dataframe)[-c(1,ncol(snps_dataframe))]){
        sub = NULL
        sub2 = NULL
        regression = NULL
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    
        if (gene_type =='v_gene'){
            regression = glm(formula = v_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_v_gene_count)
        } else if (gene_type =='d0_gene'){
            regression = glm(formula = d0_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_d_gene_count)
        } else if (gene_type =='d1_gene'){
            regression = glm(formula = d1_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_d_gene_count)
        } else if (gene_type =='j_gene'){
            regression = glm(formula = j_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_j_gene_count)
        }
        simple_regression_results = rbind(simple_regression_results, data.table(snp = snpID, intercept = as.numeric(coef(regression)[1]), slope = as.numeric(coef(regression)[2])))
        se = bootstrap_cluster(sub2[snp != "NA"], repetitions, gene_type)
        bootstrap_results = rbind(bootstrap_results, data.table(snp = snpID, standard_error = se))
    }
    together = bootstrap_regression_combine(bootstrap_results, simple_regression_results)
    return(together)
}





