args = commandArgs(trailingOnly=TRUE)
TRIM_TYPE <<- args[1]
PCA_STRUCTURE_CORRECTION = args[2]
OUTPUT_PATH = args[3]
PROJECT_PATH = args[4]

source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/config.R"))
source(paste0(PROJECT_PATH, '/tcr-gwas/trimming_regression/scripts/compile_data_functions.R'))


make_regression_file_path(TRIM_TYPE, PCA_STRUCTURE_CORRECTION)
