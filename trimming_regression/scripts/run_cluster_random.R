# import functions
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_data_functions.R")

# restrict threads
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

# This script takes four arguments (1: starting snp position, 2: trim type, 3:
# number of cpus, 4: number of bootstrap repetitions) 
args = commandArgs(trailingOnly=TRUE)

# set snp count per job
count = 1000

# set data file path
data_file_path = '/home/mrussel2/tcr-gwas/_ignore/'

# Read in snp and genotype files
snp_data = snp_file_by_snp_start(snp_start = as.numeric(args[1]), count)
genotype_data = compile_all_genotypes(snp_start = as.numeric(args[1]), count)
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)

pvalue_boot_threshold = 5e-5
pca_structure_correction = 'True'

# Run regression/bootstrap
run_snps_trimming_snp_list_cluster(snp_list = snp_data, genotype_list = genotype_data_filtered, trim_type = args[2], pca_structure_correction, pvalue_boot_threshold, write_table = 'False', ncpus = as.numeric(args[3]), maf_cutoff = 0.05, data_file_path)


print(paste0("Finished regressions for ", args[2], " for snps ", args[1], '-', as.character(as.numeric(args[1])+count)))

