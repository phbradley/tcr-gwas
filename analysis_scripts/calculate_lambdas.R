library(data.table)
library(plyr)
source('config/config.R')
source('config/file_paths.R')

args = commandArgs(trailingOnly=TRUE)

PHENOTYPE <<- args[1]
NCPU <<- as.numeric(args[2])
REPETITIONS <<- 100
CHAIN <<- 'beta'

source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))

gwas_data = compile_manhattan_plot_data(PHENOTYPE)
lambdas = get_lambdas(gwas_data, PHENOTYPE) 
