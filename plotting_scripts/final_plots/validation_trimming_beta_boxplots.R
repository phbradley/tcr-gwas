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

v_trim_data = compile_data(phenotype = 'v_trim')
j_trim_data = compile_data(phenotype = 'j_trim')

# cdr3 = combine_genes_by_common_cdr3()
# j_gene = 'TRBJ1-1*01'
# j_gene_cdr3_id = cdr3[id == j_gene]$cdr3_gene_group
# j_trim_data_subset = j_trim_data[cdr3_gene_group == j_gene_cdr3_id]

# # v_gene = 'TRBV12-3*01'
# v_gene = 'TRBV23-1*01'
# v_gene_cdr3_id = cdr3[id == v_gene]$cdr3_gene_group
# v_trim_data_subset = v_trim_data[cdr3_gene_group == v_gene_cdr3_id]

v_trim_data_mean = v_trim_data[, sum(v_trim*v_gene_count)/sum(v_gene_count), by = .(localID, productive, rs12768894, rs3762093)]
setnames(v_trim_data_mean, 'V1', 'v_trim')

j_trim_data_mean = j_trim_data[, sum(j_trim*j_gene_count)/sum(j_gene_count), by = .(localID, productive, rs12768894, rs3762093)]
setnames(j_trim_data_mean, 'V1', 'j_trim')

# Make figures for trimming: 
j_trim_boxplot = boxplot_by_snp(data = j_trim_data_mean, snp = 'rs12768894', feature_of_interest = 'j_trim')
final_j_trim_boxplot = j_trim_boxplot + xlab('rs12768894 SNP genotype') + ylab('Number of J-gene deletions') + ggtitle('')

v_trim_boxplot = boxplot_by_snp(data = v_trim_data_mean, snp = 'rs12768894', feature_of_interest = 'v_trim')
final_v_trim_boxplot = v_trim_boxplot + xlab('rs12768894 SNP genotype') + ylab('Number of V-gene deletions') + ggtitle('')

# combine v and j boxplots for alpha and beta chains

aligned_plots = align_plots(final_v_trim_boxplot, final_j_trim_boxplot, align = 'v')

together = plot_grid(aligned_plots[[1]], aligned_plots[[2]], nrow = 2, label_size = 40, labels = "AUTO", scale = 0.97)
ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/multipanel_trimming_boxplot_', CHAIN, '_rs12768894.pdf'), plot = together, width = 12, height = 17, units = 'in', dpi = 750, device = cairo_pdf)


