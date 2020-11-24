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
BY_RACE <<- args[8]

# set config variables for regression
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/config.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster_by_race.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))

dntt = find_snp_start_by_position(chromosome = 10, position1 = 98064085, position2 =98098321)
gene = get(gene_name)

snp_start = gene[1]
count = gene[2]


# Read in snp and genotype files
snp_data = snp_file_by_snp_start(snp_start, count)
genotype_data = compile_all_genotypes(snp_start, count)
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)

# Run regression/bootstrap
together_by_race = run_snps_trimming_snp_list_cluster_by_race(snp_list = snp_data, 
                                   genotype_list = genotype_data_filtered, 
                                   trim_type = TRIM_TYPE, 
                                   pca_structure_correction = 'False', 
                                   pca_type = 'none', 
                                   write_table = 'False') 
together_by_race$pca_correction = 'none'

together_all = run_snps_trimming_snp_list_cluster(snp_list = snp_data, 
                                   genotype_list = genotype_data_filtered, 
                                   trim_type = TRIM_TYPE, 
                                   pca_structure_correction = 'True', 
                                   pca_type = '8_pc_air', 
                                   write_table = 'False') 

together_all$race = 'all_together'
together_all$pca_correction = '8_pc_air'

together = rbind(together_by_race, together_all)

file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/by_race/', gene_name, '_', TRIM_TYPE, '_by_race.tsv')

write.table(as.data.frame(together), file= file_name, quote=FALSE, sep='\t', col.names= NA)

