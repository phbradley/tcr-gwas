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
source(paste0(PROJECT_PATH, '/tcr-gwas/config/validation_file_paths_', CHAIN, '.R'))

source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/validation_cohort_analysis/regression_functions_validation_cohort.R'))

genotypes = compile_all_genotypes(VALIDATION_SNP_PATH) 
phenotypes = compile_phenotype_data() 

BOOTSTRAP_PVALUE_CUTOFF <<- 1

validate = execute_regressions(genotypes, phenotypes, write.table = FALSE)
validate$one_side_pvalue = 0.5 * validate$pvalue

if (CHAIN == 'beta'){
    source('config/file_paths.R')
    source(paste0(PROJECT_PATH, '/tcr-gwas/scripts/regression_functions.R'))
    source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
    source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))

    original_data = compile_manhattan_plot_data(PHENOTYPE)
    original_data = combine_rsids(original_data)
    original_subset = original_data[substring(rsid, 1, 10) == 'rs12768894' | rsid == 'rs3762093']

    genotypes = compile_all_genotypes_snp_list(unique(original_subset$snp))
    phenotypes = compile_phenotype_data() 
    snp_meta_data = snp_file_by_snp_list(original_subset$snp)

    BOOTSTRAP_PVALUE_CUTOFF <<- 1

    original = execute_regressions(snp_meta_data, genotypes, phenotypes, write.table = FALSE)

    print(original)
} else {
    print('no original analysis for alpha chain')
}
print(validate)
