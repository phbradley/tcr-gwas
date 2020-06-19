library(data.table)
library(lme4)
library(boot)
library(ClusterBootstrap)

subset_data_snp <- function(snpID, snps_dataframe, condensed_trimming_dataframe, productive){
    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
    sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    return(sub2[snp != "NA"])
}

#model_coef <- function(data, index){
#    coef(glm(formula = v_trim ~ snp, data = data, weights = weighted_v_gene_count, subset = index))
#}
#
#bootstrap_coef <- function(data, repetitions){
#    boot = boot(data, model_coef, repetitions, weights = data$weighted_v_gene_count)
#    standard_error = sd(boot$t[,2])
#    return(standard_error)
#}

model_se <- function(data, index, gene_type){
    if (gene_type =='v_gene'){
        object = glm(formula = v_trim ~ snp, data = data, weights = weighted_v_gene_count, subset = index)
    } else if (gene_type =='d0_gene'){
        object = glm(formula = d0_trim ~ snp, data = data, weights = weighted_d_gene_count, subset = index)
    } else if (gene_type =='d1_gene'){
        object = glm(formula = d1_trim ~ snp, data = data, weights = weighted_d_gene_count, subset = index)
    } else if (gene_type =='j_gene'){
        object = glm(formula = j_trim ~ snp, data = data, weights = weighted_j_gene_count, subset = index)
    }
    coefs = summary(object)$coefficients
    coefs[,'Std. Error']
}

bootstrap_se <- function(data, repetitions, gene_type){
    if (gene_type =='v_gene'){
        weight = data$weighted_v_gene_count
    } else if (gene_type =='d0_gene' | gene_type =='d1_gene'){
        weight = data$weighted_d_gene_count
    } else if (gene_type =='j_gene'){
        weight = data$weighted_j_gene_count
    }
    boot = boot(data, model_se, repetitions, gene_type = gene_type, weights = weight, parallel="multicore", ncpus = 20)
    standard_error = colMeans(boot$t)[2]
    return(standard_error)
}

bootstrap_cluster <- function(data, repetitions, trim_type){
    if (trim_type =='v_trim'){
        boot = clusbootglm(v_trim ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='d0_trim'){
        boot = clusbootglm(d0_trim ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='d1_trim'){
        boot = clusbootglm(d1_trim ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='j_trim'){
        boot = clusbootglm(j_trim ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='vj_insert'){
        boot = clusbootglm(vj_insert ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='dj_insert'){
        boot = clusbootglm(dj_insert ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='vd_insert'){
        boot = clusbootglm(vd_insert ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    }
    standard_error = boot$boot.sds[2]
    return(standard_error)
}

regression_weighted_bootstrap_se <- function(snps_dataframe, condensed_trimming_dataframe, productive, repetitions, gene_type){
    bootstrap_results = data.frame()
    for (snpID in names(snps_dataframe)[-c(1,ncol(snps_dataframe))]){
        data = subset_data_snp(snpID, snps_dataframe, condensed_trimming_dataframe, productive)
        se = bootstrap_cluster(data, repetitions, gene_type)
        bootstrap_results = rbind(bootstrap_results, data.frame(snp = snpID, standard_error = se))
    }
    return(bootstrap_results)
}

bootstrap_regression_combine <- function(bootstrap_results, regression_results){
    if (nrow(bootstrap_results) != 0 & nrow(regression_results) != 0){
        together = merge(bootstrap_results, regression_results, by = "snp")
        together$zscore = together$slope/together$standard_error
        together$pvalue = 2*pnorm(-abs(together$zscore))
    } else {
        together = data.table()
    }
    return(together)
}


