

execute_regression <- function(snp, snps_gds, snp_id_list, trim_type, gene_type, condensing, trimming_data, repetitions, weighting, gene_conditioning, random_effects, regression_dataframe){
    # compile genotypes for snp
    snp_genotypes = compile_genotype_data(snps_gds, snp)
    index = match(snp, snp_id_list)
    # skip iteration if the genotypes are all the same...
    if (length(unique(snp_genotypes$snp)[!is.na(unique(snp_genotypes$snp))]) <= 1){
        regression_dataframe = rbind(regression_dataframe, data.frame(snp = paste0('snp',snp), intercept = 'NA', slope = 'NA', standard_error = 'NA', pvalue = 'NA', productivity = 'NA'))
        print(paste0("no regression needed for snp data for ", index, " of ", length(snp_id_list), " for ", trim_type))
    } else {
        # do regression, bootstrap
        if (condensing == "by_gene" | condensing == 'gene_cross'){
            if (random_effects == 'True'){
                regression_productive = suppressMessages(trimming_regression(snps_dataframe = snp_genotypes, condensed_trimming_dataframe = trimming_data, productive = "True",trim_type = trim_type, gene_type = gene_type, bootstrap_repetitions = repetitions,gene_conditioning,  weighting))
                regression_NOT_productive = suppressMessages(trimming_regression(snps_dataframe = snp_genotypes, condensed_trimming_dataframe = trimming_data, productive = "False", trim_type = trim_type, gene_type = gene_type, bootstrap_repetitions = repetitions,  gene_conditioning, weighting))
            } else {
                regression_productive = simple_trimming_snp_regression(snp_genotypes, trimming_data, productive = "True", trim_type =trim_type, gene_type = gene_type, weighting, gene_conditioning, python_test = 'False')
                regression_NOT_productive = simple_trimming_snp_regression(snp_genotypes, trimming_data, productive = "False", trim_type = trim_type, gene_type = gene_type, weighting, gene_conditioning, python_test = 'False')
            }  
        } else if (condensing == 'by_patient'| condensing == 'phil'){
            regression_productive = simple_trimming_snp_regression(snp_genotypes, trimming_data, productive = "True", trim_type =trim_type, weighting, gene_conditioning, python_test = 'True')
            regression_NOT_productive = simple_trimming_snp_regression(snp_genotypes, trimming_data, productive = "False", trim_type = trim_type, weighting, gene_conditioning, python_test = 'True')
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
        print(paste0("finished regression for snp data for ", index, " of ", length(snp_id_list), " for ", trim_type))
    }
    return(regression_dataframe)
}
