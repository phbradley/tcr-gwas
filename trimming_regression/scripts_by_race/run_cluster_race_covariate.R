## restrict threads
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

# This script takes four arguments (1: starting snp position, 2: trim type, 3:
# pca boolean, 4: number of cpus, 5: project_path, 6: output_path, 7: pca type
# (want `8_pca_air`)) 
args = commandArgs(trailingOnly=TRUE)

gene_name = args[1]
stopifnot(gene_name %in% c('dntt'))

TRIM_TYPE <<- args[2]
NCPU <<- args[3]
PROJECT_PATH <<- args[4]
OUTPUT_PATH <<- args[5]
BY_RACE <<- args[6]

# set config variables for regression
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts_by_race/config.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts_by_race/run_bootstrap_regression_all_snps_functions_cluster_race_covariate.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts_by_race/compile_data_functions.R"))

dntt = find_snp_start_by_position(chromosome = 10, position1 = 98064085, position2 =98098321)
gene = get(gene_name)

snp_start = gene[1]
count = gene[2]


# Read in snp and genotype files
snp_data = snp_file_by_snp_start(snp_start, count)
genotype_data = compile_all_genotypes(snp_start, count)
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)

together = run_snps_trimming_snp_list_cluster_race_covariate(snp_list = snp_data, 
                                   genotype_list = genotype_data_filtered, 
                                   trim_type = TRIM_TYPE, 
                                   pca_structure_correction = 'False', 
                                   pca_type = 'none', 
                                   by_race = BY_RACE,
                                   write_table = 'False') 

file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/by_race/', gene_name, '_', TRIM_TYPE, '_with_race_covariate.tsv')

write.table(as.data.frame(together), file= file_name, quote=FALSE, sep='\t', col.names= NA)

