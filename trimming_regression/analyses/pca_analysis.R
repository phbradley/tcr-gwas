library(data.table)
setDTthreads(1)


source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_regression_data_functions.R')

ethnicity = fread(file = '/home/mrussel2/tcr-gwas/_ignore/snp_data/ethnicity_data.csv') 

pc = read_genotype_pca()


