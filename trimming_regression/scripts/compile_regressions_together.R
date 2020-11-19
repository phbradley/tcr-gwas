args = commandArgs(trailingOnly=TRUE)

TRIM_TYPE = args[1]
PCA_STRUCTURE_CORRECTION = args[2]
PCA_TYPE = args[3]
PROJECT_PATH = args[4]
OUTPUT_PATH = args[5]

source(paste0(PROJECT_PATH, '/tcr-gwas/trimming_regression/scripts/compile_regression_functions.R'))

# Use this, for running different trimming type compiles in parallel, but each
# is indidividually run sequentially

compile_all_data_from_cluster_sequential(TRIM_TYPE, PCA_STRUCTURE_CORRECTION, PCA_TYPE) 
