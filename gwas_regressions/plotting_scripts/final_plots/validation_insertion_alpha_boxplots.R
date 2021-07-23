library(cowplot)
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(Cairo)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)

CHAIN <<- 'alpha'
source('config/config.R')
source(paste0('config/validation_file_paths_', CHAIN, '.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/validation_plot_functions.R'))

vj_insertion_data = compile_data(phenotype = 'vj_insert')


# Make figures for insertions: 
vj_insert_boxplot = boxplot_by_snp(data = vj_insertion_data, snp = 'rs3762093', feature_of_interest = 'vj_insert')
final_vj_insert_boxplot = vj_insert_boxplot + xlab('rs3762093 SNP genotype') + ylab('Number of V-J-gene junction N-inserts') + ggtitle('')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/vj_insert_boxplot_', CHAIN, '_rs3762093.pdf'), plot = final_vj_insert_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)



