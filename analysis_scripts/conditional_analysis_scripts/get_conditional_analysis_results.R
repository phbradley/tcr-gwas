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

source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
original_analysis = compile_manhattan_plot_data(c(PHENOTYPE)) 
gene = GENE_ANNOTATIONS[gene_common_name == GENE]

original_analysis_gene = original_analysis[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
significance_cutoff = determine_significance_cutoff(0.05, type = GENE, original_analysis_gene)

original_analysis_gene_sig = original_analysis_gene[pvalue < significance_cutoff][order(pvalue)]

source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/conditional_analysis_scripts/conditional_regression_functions.R'))


files = fs::dir_ls(path = make_regression_file_path())
together = data.table()
for (file in files){
    together = rbind(together, fread(file), fill=TRUE)
}

results = together[pvalue < significance_cutoff]
print("conditional analysis results for productive and non-productive cases:")
print(results)

print("original analysis output for conditional ananlysis sig. snp:")
if (nrow(results) != 0){
    print(original_analysis_gene_sig[snp %in% unique(results$snp)])
}

print("original analysis output for most significant snp:")
print(original_analysis_gene_sig[snp %in% unique(together$conditional_snps)])
