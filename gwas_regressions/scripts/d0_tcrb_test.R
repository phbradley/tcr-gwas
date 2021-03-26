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


START <<- 16730300
most_sig_snp = 16814536
PHENOTYPE <<- 'd0_trim'
NCPU <<- 4

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))

SNPS_PER_JOB =100
genotypes = compile_all_genotypes(as.numeric(START), as.numeric(SNPS_PER_JOB))
phenotypes = compile_phenotype_data() 
snp_meta_data = snp_file_by_snp_start(as.numeric(START), as.numeric(SNPS_PER_JOB))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/phenotype_classes/snp_interaction_trimming_class_functions.R'))
with_interaction_test = execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = FALSE)

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/phenotype_classes/trimming_class_functions.R'))
normal_test = execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = FALSE)


#Now test artemis snps for j_trim in same scheme

START <<- 20615300
most_sig_snp = 20717772
PHENOTYPE <<- 'j_trim'
NCPU <<- 4

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))

SNPS_PER_JOB =100
genotypes = compile_all_genotypes(as.numeric(START), as.numeric(SNPS_PER_JOB))
phenotypes = compile_phenotype_data() 
snp_meta_data = snp_file_by_snp_start(as.numeric(START), as.numeric(SNPS_PER_JOB))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/phenotype_classes/snp_interaction_trimming_class_functions.R'))
with_interaction_test_artemis = execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = FALSE)

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/phenotype_classes/trimming_class_functions.R'))
normal_test_artemis = execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = FALSE)

