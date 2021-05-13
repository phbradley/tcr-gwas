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

args = commandArgs(trailingOnly=TRUE)

PHENOTYPE <<- args[1]
NCPU <<- as.numeric(args[2])
CHAIN <<- args[3]

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/config/validation_file_paths_', CHAIN, '.R'))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/validation_cohort_analysis/regression_functions_validation_cohort.R'))

genotypes = compile_all_genotypes(VALIDATION_SNP_PATH) 
phenotypes = compile_phenotype_data() 

execute_regressions(genotypes, phenotypes, write.table = TRUE)

