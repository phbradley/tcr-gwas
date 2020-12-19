source('config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/file_paths.R'))

library(data.table)
setDTthreads(1)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(foreach)
library(doParallel)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

PHENOTYPE <<- args[1]
NCPU <<- args[2]

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/src/regression_functions.R'))

compile_condensed_tcr_repertoire_data()
