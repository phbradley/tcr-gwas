library(lme4)
library(data.table)

source("trimming_bootstrap_functions.R")


trimming_snp_regression_weighted_varying_int_subject <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, repetitions){
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
    bonferroni = 0.05/35481497
    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    simple_regression_results = data.table()
    bootstrap_results = data.table()
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    for (snpID in names(snps_dataframe)[-c(1)]){
        sub = NULL
        sub2 = NULL
        regression = NULL
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    
        if (trim_type =='v_trim'){
            regression = lmer(formula = v_trim ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_v_gene_count, control=control)
        } else if (trim_type =='d0_trim'){
            regression = lmer(formula = d0_trim ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_d_gene_count, control=control)
        } else if (trim_type =='d1_trim'){
            regression = lmer(formula = d1_trim ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_d_gene_count, control=control)
        } else if (trim_type =='j_trim'){
            regression = lmer(formula = j_trim ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_j_gene_count, control=control)
        } else if (trim_type =='vj_insert'){
            regression = lmer(formula = vj_insert ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_vj_gene_count, control=control)
        } else if (trim_type =='dj_insert'){
            regression = lmer(formula = dj_insert ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_dj_gene_count, control=control)
        } else if (trim_type =='vd_insert'){
            regression = lmer(formula = vd_insert ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_vd_gene_count, control=control)
        }
        simple_regression_results = rbind(simple_regression_results, data.table(snp = snpID, intercept = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][1])))), slope = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][2]))))))
        
        se = clusboot_lmer(regression, data = sub2[snp != "NA"], cluster_variable = sub2[snp != "NA"]$localID, repetitions)[2,2]
        
        bootstrap_results = rbind(bootstrap_results, data.table(snp = snpID, standard_error = se))
    }
    together = bootstrap_regression_combine(bootstrap_results, simple_regression_results, bonferroni)
    return(together)
}



trimming_snp_regression_weighted <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, repetitions){
    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    simple_regression_results = data.table()
    bootstrap_results = data.table()
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    for (snpID in names(snps_dataframe)[-c(1)]){
        sub = NULL
        sub2 = NULL
        regression = NULL
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    
        if (trim_type =='v_trim'){
            regression = glm(formula = v_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_v_gene_count)
        } else if (trim_type =='d0_trim'){
            regression = glm(formula = d0_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_d_gene_count)
        } else if (trim_type =='d1_trim'){
            regression = glm(formula = d1_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_d_gene_count)
        } else if (trim_type =='j_trim'){
            regression = glm(formula = j_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_j_gene_count)
        } else if (trim_type =='vj_insert'){
            regression = glm(formula = vj_insert ~ snp, data = sub2[snp != "NA"], weights = weighted_vj_gene_count)
        } else if (trim_type =='dj_insert'){
            regression = glm(formula = dj_insert ~ snp, data = sub2[snp != "NA"], weights = weighted_dj_gene_count)
        } else if (trim_type =='vd_insert'){
            regression = glm(formula = vd_insert ~ snp, data = sub2[snp != "NA"], weights = weighted_vd_gene_count)
        }
        simple_regression_results = rbind(simple_regression_results, data.table(snp = snpID, intercept = as.numeric(coef(regression)[1]), slope = as.numeric(coef(regression)[2])))
        se = bootstrap_cluster(sub2[snp != "NA"], repetitions, trim_type)
        bootstrap_results = rbind(bootstrap_results, data.table(snp = snpID, standard_error = se))
    }
    together = bootstrap_regression_combine(bootstrap_results, simple_regression_results)
    return(together)
}





