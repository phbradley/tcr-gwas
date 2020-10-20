
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/make_snp_file.R")

library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)


count = 1000
# Read in snp list
snp_data = snp_file_by_snp_start(snp_start = as.numeric(args[1]), count)
genotype_data = compile_all_genotypes(snp_start = as.numeric(args[1]), count)
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)

# Run regression/bootstrap
run_snps_trimming_snp_list_cluster(snp_list = snp_data, genotype_list = genotype_data_filtered, trim_type = args[2], gene_type = 'same', condensing = 'by_gene', gene_conditioning = 'True', weighting = 'True', random_effects = 'True', repetitions = as.numeric(args[4]), pca_structure_correction = 'True', write_table = 'True', ncpus = as.numeric(args[3]), d_infer = 'True', maf_cutoff = 0.05)
print(paste0("Finished regressions for ", args[2], " for snps ", args[1], '-', as.character(as.numeric(args[1])+count)))

