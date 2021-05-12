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

START <<- 20615500 
PHENOTYPE <<- 'v_trim' 
NCPU <<- as.numeric(5)

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))

genotypes = compile_all_genotypes(as.numeric(START), as.numeric(SNPS_PER_JOB))
phenotypes = compile_phenotype_data() 
snp_meta_data = snp_file_by_snp_start(as.numeric(START), as.numeric(SNPS_PER_JOB))

results = execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = FALSE)


KEEP_MISSING_D_GENE = 'False'

no_missing_d_phenotypes = compile_phenotype_data()

no_missing_d_results = execute_regressions(snp_meta_data, genotypes, no_missing_d_phenotypes, write.table = FALSE)


pdf(file="test_d_infer_comparison_v_trim.pdf")
plot(-log10(no_missing_d_results$pvalue), -log10(results$pvalue), main = 'V-gene trimming p-value comparison\nincluding missing D-gene observations versus excluding them', xlab = '-log10(p-value) when exluding missing D-gene observations', ylab = '-log10(p-value) when including D-gene observations', ylim = c(0, 35), xlim = c(0,35))
abline(a = 0, b = 1)
dev.off()


START <<- 16730000 
PHENOTYPE <<- 'd0_trim' 
KEEP_MISSING_D_GENE = 'True'

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))

genotypes = compile_all_genotypes(as.numeric(START), as.numeric(SNPS_PER_JOB))
phenotypes = compile_phenotype_data() 
snp_meta_data = snp_file_by_snp_start(as.numeric(START), as.numeric(SNPS_PER_JOB))

results = execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = FALSE)


KEEP_MISSING_D_GENE = 'False'

no_missing_d_phenotypes = compile_phenotype_data()

no_missing_d_results = execute_regressions(snp_meta_data, genotypes, no_missing_d_phenotypes, write.table = FALSE)


pdf(file="test_d_infer_comparison_d0_trim.pdf")
plot(-log10(no_missing_d_results$pvalue), -log10(results$pvalue), main = 'D0-gene trimming p-value comparison\ninferring missing D-gene observations versus excluding them', xlab = '-log10(p-value) when exluding missing D-gene observations', ylab = '-log10(p-value) when inferring D-gene observations')
abline(a = 0, b = 1)
dev.off()




