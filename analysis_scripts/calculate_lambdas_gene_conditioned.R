source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/config/file_paths.R'))

library(data.table)
setDTthreads(1)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(foreach)
library(doParallel)
library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

PHENOTYPE <<- args[1]
NCPU <<- as.numeric(args[2])
REPETITIONS <<- 100
CHAIN <<- 'beta'

source(paste0(PROJECT_PATH, '/tcr-gwas/scripts/regression_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

results = get_random_boostrap_results()
lambdas = get_lambdas(results, PHENOTYPE)
