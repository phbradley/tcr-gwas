
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

execute_regression <- function(snp, snp_list, genotype_list, trim_type, gene_type, condensing, trimming_data, repetitions, weighting, gene_conditioning, random_effects, regression_dataframe){
    genotypes_temp = as.data.frame(genotype_list[, as.character(snp)])
    snp_genotypes = data.frame(rownames(genotypes_temp),genotypes_temp)
    colnames(snp_genotypes) = c('localID', paste0('snp',snp))
    rownames(snp_genotypes) = NULL

    index = which(as.numeric(colnames(genotype_list)) == snp)

    # skip iteration if the genotypes are all the same...
    if (length(unique(snp_genotypes$snp)[!is.na(unique(snp_genotypes$snp))]) <= 1){
        temp_regression_dataframe = merge(snp_list, data.frame(snp = snp, intercept = 'NA', slope = 'NA', standard_error = 'NA', pvalue = 'NA', productivity = 'NA'), by = 'snp')
        temp_regression_dataframe$snp = paste0('snp', temp_regression_dataframe$snp)
        regression_dataframe = rbind(regression_dataframe, temp_regression_dataframe)
        print(paste0("no regression needed for snp data for ", index, " of ", ncol(genotype_list), " for ", trim_type))
    } else {
        # do regression, bootstrap
        if (condensing == "by_gene" | condensing == 'gene_cross'){
            if (random_effects == 'True'){
                regression_productive = suppressMessages(trimming_regression(snps_dataframe = snp_genotypes, condensed_trimming_dataframe = trimming_data, productive = "True",trim_type = trim_type, gene_type = gene_type, bootstrap_repetitions = repetitions,gene_conditioning,  weighting, snp_list))
                regression_NOT_productive = suppressMessages(trimming_regression(snps_dataframe = snp_genotypes, condensed_trimming_dataframe = trimming_data, productive = "False", trim_type = trim_type, gene_type = gene_type, bootstrap_repetitions = repetitions,  gene_conditioning, weighting, snp_list))
            } else {
                regression_productive = simple_trimming_snp_regression(snp_genotypes, trimming_data, productive = "True", trim_type =trim_type, gene_type = gene_type, weighting, gene_conditioning, python_test = 'False', snp_list)
                regression_NOT_productive = simple_trimming_snp_regression(snp_genotypes, trimming_data, productive = "False", trim_type = trim_type, gene_type = gene_type, weighting, gene_conditioning, python_test = 'False', snp_list)
            }  
        } else if (condensing == 'by_patient'| condensing == 'phil'){
            regression_productive = simple_trimming_snp_regression(snp_genotypes, trimming_data, productive = "True", trim_type =trim_type, weighting, gene_conditioning, python_test = 'True', snp_list)
            regression_NOT_productive = simple_trimming_snp_regression(snp_genotypes, trimming_data, productive = "False", trim_type = trim_type, weighting, gene_conditioning, python_test = 'True', snp_list)
        } 

        # if regression, bootstrap conducted, add them to the results tables
        if (nrow(regression_productive) != 0){
            regression_productive$productivity = 'productive'
            regression_dataframe = rbind(regression_dataframe, regression_productive)
        }
        if (nrow(regression_NOT_productive) != 0){
            regression_NOT_productive$productivity = 'NOT_productive'
            regression_dataframe = rbind(regression_dataframe, regression_NOT_productive)
        }
        print(paste0("finished regression for snp data for ", index, " of ", ncol(genotype_list), " for ", trim_type))
    }
    return(regression_dataframe)
}