# This script will combine all genome-wide regression runs (1000 per file) into one master file
source('config/config.R')

library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(vroom)
library(fs)
library(data.table)
setDTthreads(1)

args = commandArgs(trailingOnly=TRUE)

PHENOTYPE <<- args[1]
NCPU <<- as.numeric(args[2])

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))

compile_regressions()
