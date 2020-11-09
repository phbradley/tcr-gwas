library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

# This function actually executes the regressions on the cluster

execute_regression <- function(snp, list_of_snps, snp_list, genotype_list, trim_type, condensing, trimming_data, repetitions, pca_structure_correction, weighting, gene_conditioning, random_effects, regression_dataframe){
    if (as.character(snp) %in% colnames(genotype_list)){
        genotypes_temp = as.data.frame(genotype_list[, as.character(snp)])
        snp_genotypes = data.frame(rownames(genotype_list),genotypes_temp)
        colnames(snp_genotypes) = c('localID', paste0('snp',snp))
        rownames(snp_genotypes) = NULL

        index = which(as.numeric(list_of_snps) == snp)

        # skip iteration if the genotypes are all the same...
        if (length(unique(snp_genotypes$snp)[!is.na(unique(snp_genotypes$snp))]) <= 1){
            temp_regression_dataframe = merge(snp_list, data.frame(snp = snp, intercept = 'NA', slope = 'NA', standard_error = 'NA', pvalue = 'NA', productivity = 'NA'), by = 'snp')
            temp_regression_dataframe$snp = paste0('snp', temp_regression_dataframe$snp)
            regression_dataframe = rbind(regression_dataframe, temp_regression_dataframe)
            print(paste0("no regression needed for snp data for ", index, " of ", ncol(genotype_list), " for ", trim_type))
        } else {
            # do regression, bootstrap
            regression_productive = suppressMessages(trimming_regression(snps_dataframe = snp_genotypes, condensed_trimming_dataframe = trimming_data, condensing, productive = "True", trim_type = trim_type, bootstrap_repetitions = repetitions, pca_structure_correction, gene_conditioning, weighting, random_effects, snp_list))
            regression_NOT_productive = suppressMessages(trimming_regression(snps_dataframe = snp_genotypes, condensed_trimming_dataframe = trimming_data,condensing, productive = "False", trim_type = trim_type, bootstrap_repetitions = repetitions,pca_structure_correction, gene_conditioning, weighting, random_effects, snp_list))
            
            if (nrow(regression_productive) != 0){
                regression_productive$productivity = 'productive'
                regression_dataframe = rbind(regression_dataframe, regression_productive)
            }
            if (nrow(regression_NOT_productive) != 0){
                regression_NOT_productive$productivity = 'NOT_productive'
                regression_dataframe = rbind(regression_dataframe, regression_NOT_productive)
            }
            print(paste0("finished regression for snp data for ", index, " of ", length(list_of_snps), " for ", trim_type))
        }
    }
    return(regression_dataframe)
}


