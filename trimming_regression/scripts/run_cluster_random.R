# import functions
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/make_snp_file.R")

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

# parse type (either 'insert' or 'trim')
stopifnot(args[2] %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'vj_insert'))
type = strsplit(args[2], '_')[[1]][2]


# set variables conditional on type (i.e. if we are doing insertion
# regressions, we want to condense the trimming data by patient without
# inferring d_gene, and regress not conditioning out TCRB gene and without 
# random effects)
condensing = ifelse(type == 'insert', 'by_patient', 'by_gene')
gene_conditioning = ifelse(type == 'insert', 'False', 'True')
random_effects = ifelse(type == 'insert', 'False', 'True')
d_infer = ifelse(type == 'insert', 'False', 'True')

condensing = ifelse(type == 'insert', 'by_patient', 'by_patient')
gene_conditioning = ifelse(type == 'insert', 'False', 'False')
random_effects = ifelse(type == 'insert', 'False', 'False')
d_infer = ifelse(type == 'insert', 'False', 'False')
# set pca_structure_correction variable (only want pca structure correction for
# second pass analysis (i.e. repetitions > 0)
pca_structure_correction = ifelse(as.numeric(args[4]) == 0, 'False', 'True')


# Run regression/bootstrap
run_snps_trimming_snp_list_cluster(snp_list = snp_data, genotype_list = genotype_data_filtered, trim_type = args[2], gene_type = 'same', condensing, gene_conditioning, weighting = 'True', random_effects, repetitions = as.numeric(args[4]), pca_structure_correction, write_table = 'True', ncpus = as.numeric(args[3]), d_infer, maf_cutoff = 0.05, data_file_path)


print(paste0("Finished regressions for ", args[2], " for snps ", args[1], '-', as.character(as.numeric(args[1])+count)))

