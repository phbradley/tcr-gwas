library(data.table)
library(lme4)
library(boot)

subset_data_snp <- function(snpID, snps_dataframe, condensed_trimming_dataframe, productive){
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    sub = data.frame(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
    sub2 = as.data.table(merge(sub, condensed_trimming_dataframe, by = "localID"))
    return(sub2[snp != "NA"])
}

model_coef <- function(data, index){
    coef(glm(formula = v_trim ~ snp, data = data, weights = weighted_v_gene_count, subset = index))
}

bootstrap_coef <- function(data, repetitions){
    boot = boot(data, model_coef, repetitions, weights = data$weighted_v_gene_count)
    standard_error = sd(boot$t[,2])
    return(standard_error)
}

model_se <- function(data, index){
    object = glm(formula = v_trim ~ snp, data = data, weights = weighted_v_gene_count, subset = index)
    coefs = summary(object)$coefficients
    coefs[,'Std. Error']
}

bootstrap_se <- function(data, repetitions){
    boot = boot(data, model_se, repetitions, weights = data$weighted_v_gene_count)
    standard_error = colMeans(boot$t)[2]
    return(standard_error)
}

regression_weighted_bootstrap_se <- function(snps_dataframe, condensed_trimming_dataframe, productive, repetitions){
    bootstrap_results = data.frame()
    for (snpID in names(snps_dataframe)[-c(1,1185)]){
        data = subset_data_snp(snpID, snps_dataframe, condensed_trimming_dataframe, productive)
        se = bootstrap_se(data, repetitions)
        bootstrap_results = rbind(bootstrap_results, data.frame(snp = snpID, standard_error = se))
    }
    return(bootstrap_results)
}

bootstrap_regression_combine <- function(bootstrap_results, regression_results){
    together = merge(bootstrap_results, regression_results, by = "snp")
    together$zscore = together$slope/together$standard_error
    together$pvalue = dnorm(together$zscore)
    return(together)
}


