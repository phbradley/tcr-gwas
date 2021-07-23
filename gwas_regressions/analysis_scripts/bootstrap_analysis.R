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
library(plyr)

args = commandArgs(trailingOnly=TRUE)

PHENOTYPE <<- args[1]
NCPU <<- as.numeric(args[2])

CHAIN <<- 'beta'

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
original_analysis = compile_manhattan_plot_data(c(PHENOTYPE)) 
snps = unique(original_analysis$snp)

random_snps = sample(1:length(snps), 100, replace = T)

random_subset = snps[random_snps]

genotypes = compile_all_genotypes_snp_list(random_subset)
phenotypes = compile_phenotype_data() 
snp_meta_data = snp_file_by_snp_list(random_subset)

BOOTSTRAP_PVALUE_CUTOFF <<- 1
REPETITIONS <<- 100 

bootstrap_subset = execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = FALSE)

path = paste0(OUTPUT_PATH, '/results/bootstrap_lambda_analysis/', PHENOTYPE, '_', REPETITIONS)
dir.create(path)

fwrite(bootstrap_subset, file = paste0(OUTPUT_PATH, '/results/bootstrap_lambda_analysis/', PHENOTYPE, '_', REPETITIONS, '/bootstrap_subset_', PHENOTYPE, '_', REPETITIONS, '_bootstraps_', random_subset[1], '.tsv'))
