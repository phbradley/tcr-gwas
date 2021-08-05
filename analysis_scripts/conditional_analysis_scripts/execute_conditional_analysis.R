source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/config/file_paths.R'))

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

GENE <<- args[1]
PHENOTYPE <<- args[2]
NCPU <<- as.numeric(args[3])

source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
original_analysis = compile_manhattan_plot_data(c(PHENOTYPE)) 
gene = GENE_ANNOTATIONS[gene_common_name == GENE]

original_analysis_gene = original_analysis[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
# multiply by two since we are only looking at one productivity status at a time
if (PHENOTYPE == 'gene_usage'){
    significance_cutoff = determine_significance_cutoff(0.05, type = GENE, original_analysis_gene, phenotype = unique(original_analysis_gene$gene))*2
} else {
    significance_cutoff = determine_significance_cutoff(0.05, type = GENE, original_analysis_gene)*2
}
original_analysis_gene_sig = original_analysis_gene[pvalue < significance_cutoff][order(pvalue)]

source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/conditional_analysis_scripts/conditional_regression_functions.R'))


for (productivity in c(TRUE, FALSE)){
         genotypes = compile_all_genotypes_snp_list(unique(original_analysis_gene[productive == productivity]$snp))
         phenotypes = compile_phenotype_data()
         # snp_meta_data = snp_file_by_snp_list(original_analysis_gene[productive == productivity]$snp)
         
         if (nrow(original_analysis_gene_sig[productive == productivity]) > 0){
            conditional_snp_start = original_analysis_gene_sig[productive == productivity][order(pvalue)][1]$snp
            execute_conditional_regressions(genotypes, phenotypes, conditional_snp_start, productivity, significance_cutoff)
         } else {
            print(paste0('No significant snps in the ', GENE, ' region for ', PHENOTYPE, ' of ', productivity, ' productivity TCRs'))
         }
}

write_final_compiled_file(significance_cutoff)
