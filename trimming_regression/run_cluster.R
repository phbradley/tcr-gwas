
source("/home/mrussel2/tcr-gwas/trimming_regression/run_bootstrap_regression_all_snps_functions_cluster.R")
#source("find_correlated_snps.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/make_snp_file.R")

args = commandArgs(trailingOnly=TRUE)


count = 5000
# Read in snp list
snp_data = snp_file_by_snp_start(snp_start = as.numeric(args[1]), count)
genotype_data = compile_all_genotypes(snp_start = as.numeric(args[1]), count)
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)

# Run regression/bootstrap
run_snps_trimming_snp_list_cluster(snp_list = snp_data, genotype_list = genotype_data_filtered, trim_type = args[2], gene_type = 'same', condensing = 'by_gene', gene_conditioning = 'True', weighting = 'True', random_effects = 'False', repetitions = 0, write_table = 'True', ncpus = 4)
print(paste0("Finished regressions for ", args[2], " for snps ", args[1], '-', as.character(as.numeric(args[1])+count)))

