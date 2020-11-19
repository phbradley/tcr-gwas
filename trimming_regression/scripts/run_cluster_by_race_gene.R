# This script takes four arguments (1: starting snp position, 2: trim type, 3:
# number of cpus, 4: project_path, 5: output_path) 
args = commandArgs(trailingOnly=TRUE)


gene_name = args[1]
stopifnot(gene_name %in% c('dntt'))

trimming_type = args[2]
pca = args[3]
ncpu = args[4]
project_path = args[5]
output_path = args[6]
pcatype = args[7]


# import functions
source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R"))
source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster_by_race.R"))
source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))


dntt = find_snp_start_by_position(chromosome = 10, position1 = 98064085, position2 =98098321)
gene = get(gene_name)

# restrict threads
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)


# Read in snp and genotype files
snp_data = snp_file_by_snp_start(snp_start = gene[1], count = gene[2])
genotype_data = compile_all_genotypes(snp_start = gene[1], count = gene[2])
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)

# Run regression/bootstrap
together_by_race = run_snps_trimming_snp_list_cluster_by_race(snp_list = snp_data, genotype_list = genotype_data_filtered, trim_type = trimming_type, pca_structure_correction = pca, pca_type = pcatype, write_table = 'True', ncpus = as.numeric(ncpu), maf_cutoff = 0.05, random_effects)
# Note: regressions using all races are using pca population structure
# correction including Dave's 8 PCAir PCs as a controlled pvalue
together_all = run_snps_trimming_snp_list_cluster(snp_list = snp_data, genotype_list = genotype_data_filtered, trim_type = trimming_type, pca_structure_correction = 'True', pca_type = '8_pc_air', write_table = 'False', ncpus = as.numeric(ncpu), maf_cutoff = 0.05, random_effects)

together_all$race = 'all_together'

together = rbind(together_by_race, together_all)

file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/by_race/', gene_name, '_', trimming_type, '_pca-', pca, '_', pcatype, '_by_race.tsv')

write.table(as.data.frame(together), file= file_name, quote=FALSE, sep='\t', col.names= NA)

