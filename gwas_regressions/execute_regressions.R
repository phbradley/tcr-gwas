source(paste0(getwd(), '/config.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/file_paths.R'))

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
library(modules)
#TODO uncomment these
# args = commandArgs(trailingOnly=TRUE)

# START <<- args[1]
# PHENOTYPE <<- args[2]
# NCPU <<- args[3]

START = 10000000
PHENOTYPE = 'd0_trim'
NCPU = 1

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/regression_functions.R'))
genotypes = compile_all_genotypes(START, SNPS_PER_JOB)
phenotypes = compile_phenotype_data() 

execute_regressions(genotypes, phenotypes, write.table = TRUE)

print(paste0('Finished regressions for ', PHENOTYPE, ' for snps ', START, '-', as.character(as.numeric(START) + as.numeric(SNPS_PER_JOB))))
