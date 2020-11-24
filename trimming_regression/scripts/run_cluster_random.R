# restrict threads
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

# This script takes four arguments (1: starting snp position, 2: trim type, 3:
# pca boolean, 4: number of cpus, 5: project_path, 6: output_path, 7: pca type
# (want `8_pca_air`)) 
args = commandArgs(trailingOnly=TRUE)

START <<- args[1]
TRIM_TYPE <<- args[2]
PCA_STRUCTURE_CORRECTION <<- args[3]
NCPU <<- args[4]
PROJECT_PATH <<- args[5]
OUTPUT_PATH <<- args[6]
PCA_TYPE <<- args[7]

# set config variables for regression
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/config.R"))

# import functions
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))

# set snp count per job
count = 1000

# Read in snp and genotype files
snp_data = snp_file_by_snp_start(snp_start = as.numeric(START), count)
genotype_data = compile_all_genotypes(snp_start = as.numeric(START), count)
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)

# Run regression/bootstrap
run_snps_trimming_snp_list_cluster(snp_list = snp_data, 
                                   genotype_list = genotype_data_filtered, 
                                   trim_type = TRIM_TYPE, 
                                   pca_structure_correction = PCA_STRUCTURE_CORRECTION, 
                                   pca_type = PCA_TYPE, 
                                   write_table = 'True') 


print(paste0("Finished regressions for ", TRIM_TYPE, " for snps ", START, '-', as.character(as.numeric(START)+count)))

