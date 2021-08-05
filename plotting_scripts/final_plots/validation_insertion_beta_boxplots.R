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

CHAIN <<- 'beta'
source('config/config.R')
source(paste0('config/validation_file_paths_', CHAIN, '.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/validation_plot_functions.R'))

vd_insertion_data = compile_data(phenotype = 'vd_insert')
dj_insertion_data = compile_data(phenotype = 'dj_insert')


# Make figures for insertions: 
vd_insert_boxplot = boxplot_by_snp(data = vd_insertion_data, snp = 'rs3762093', feature_of_interest = 'vd_insert')
final_vd_insert_boxplot = vd_insert_boxplot + xlab('rs3762093 SNP genotype') + ylab('Number of V-D-gene junction N-inserts') + ggtitle('')

dj_insert_boxplot = boxplot_by_snp(data = dj_insertion_data, snp = 'rs3762093', feature_of_interest = 'dj_insert')
final_dj_insert_boxplot = dj_insert_boxplot + xlab('rs3762093 SNP genotype') + ylab('Number of D-J-gene junction N-inserts') + ggtitle('')

aligned_plots = align_plots(final_vd_insert_boxplot, final_dj_insert_boxplot, align = 'v')

together = plot_grid(aligned_plots[[1]], aligned_plots[[2]], nrow = 2, label_size = 40, labels = "AUTO", scale = 0.97)
ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/multipanel_insertion_boxplot_', CHAIN, '_rs3762093.pdf'), plot = together, width = 12, height = 17, units = 'in', dpi = 750, device = cairo_pdf)



