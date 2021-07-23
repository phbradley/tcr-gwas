# This script will run regressions for all snps in the indicated snp-chunk
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

args = commandArgs(trailingOnly=TRUE)

START <<- args[1]
PHENOTYPE <<- args[2]
NCPU <<- as.numeric(args[3])
CHAIN <<- 'beta'

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))

genotypes = compile_all_genotypes(as.numeric(START), as.numeric(SNPS_PER_JOB))
phenotypes = compile_phenotype_data() 
snp_meta_data = snp_file_by_snp_start(as.numeric(START), as.numeric(SNPS_PER_JOB))


execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = TRUE)

print(paste0('Finished regressions for ', PHENOTYPE, ' for snps ', START, '-', as.character(as.numeric(START) + as.numeric(SNPS_PER_JOB))))
