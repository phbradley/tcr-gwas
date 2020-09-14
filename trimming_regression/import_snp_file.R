# import snp and genotype files...they are currently in 100000 snp chunks, so need to read in the right one and extract the correct snps

get_snp_file <- function(snp_start, count){
    snp_data = read.table(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/snp_data/', snp_start, '_', count, '.txt'), sep = "\t", fill=TRUE, header = TRUE)[-1]
    return(snp_data)
}

get_genotype_file <- function(snp_start, count){
    genotype_data = as.matrix(read.table(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/genotype_data/', snp_start, '_', count, '.txt'), sep = "\t", fill=TRUE, header = TRUE, row.names=1, check.names = FALSE))
    return(genotype_data)
}