# This script will return the output regression file path
source('config/config.R')

args = commandArgs(trailingOnly=TRUE)

PHENOTYPE <<- args[1]

source(paste0(PROJECT_PATH, '/tcr-gwas/scripts/regression_functions.R'))

find_regression_file_path_for_shell()
