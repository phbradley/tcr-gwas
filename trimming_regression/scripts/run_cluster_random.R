# This script takes four arguments (1: starting snp position, 2: trim type, 3:
# number of cpus, 4: project_path, 5: output_path) 
args = commandArgs(trailingOnly=TRUE)

start = args[1]
trimming_type = args[2]
pca = args[3]
ncpu = args[4]
project_path = args[5]
output_path = args[6]

# import functions
source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R"))
source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))

# restrict threads
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

# set snp count per job
count = 1000

# Read in snp and genotype files
snp_data = snp_file_by_snp_start(snp_start = as.numeric(start), count)
genotype_data = compile_all_genotypes(snp_start = as.numeric(start), count)
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)

pvalue_boot_threshold = 5e-5

# Run regression/bootstrap
run_snps_trimming_snp_list_cluster(snp_list = snp_data, genotype_list = genotype_data_filtered, trim_type = trimming_type, pca_structure_correction = pca, pvalue_boot_threshold, write_table = 'True', ncpus = as.numeric(ncpu), maf_cutoff = 0.05)


print(paste0("Finished regressions for ", trimming_type, " for snps ", start, '-', as.character(as.numeric(start)+count)))

