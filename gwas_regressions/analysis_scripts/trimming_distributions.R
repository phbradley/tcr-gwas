source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/config/file_paths.R'))

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

PHENOTYPE <<- 'd0_trim_no_missing_d'
NCPU <<- 1 

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))

genotypes = compile_all_genotypes(as.numeric(START), as.numeric(SNPS_PER_JOB))
phenotypes = compile_phenotype_data() 

