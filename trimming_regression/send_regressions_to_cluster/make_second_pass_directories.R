
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/make_snp_file.R")

library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

all_snps = open_regressed_file_and_subset_by_pval(significance_cutoff = 5e-5, trim_type = args[1], random_effects = 'False', condensing = 'by_patient', maf_cutoff = 0.05)
count = 10

index_count = ifelse(nrow(all_snps)/count < 1, 1, nrow(all_snps)/count)
index_count = ceiling(index_count)
system(paste0('./make_combination_directories_second_pass.sh ', args[1], ' ', index_count, ' second_pass_directories'))
