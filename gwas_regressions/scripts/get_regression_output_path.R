source('config/config.R')

args = commandArgs(trailingOnly=TRUE)

PHENOTYPE <<- args[1]

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/src/regression_functions.R'))

find_regression_file_path_for_shell()
