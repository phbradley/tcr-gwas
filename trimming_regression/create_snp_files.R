source("/home/mrussel2/tcr-gwas/trimming_regression/make_snp_file.R")

count = 1000

for (snp_start in seq(01, 35500000, count)){
    snp_data = snp_file_by_snp_start(snp_start = as.numeric(snp_start), count)
    genotype_data = compile_all_genotypes(snp_start = as.numeric(snp_start), count)
    write.table(snp_data, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/snp_data/', snp_start, '_', count, '.txt'), quote=FALSE, sep='\t', col.names = NA)
    write.table(genotype_data, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/genotype_data/', snp_start, '_', count, '.txt'), quote=FALSE, sep='\t', col.names = NA)
}


