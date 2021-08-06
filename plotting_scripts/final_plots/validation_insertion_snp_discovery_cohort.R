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

args = commandArgs(trailingOnly=TRUE)

NCPU <<- as.numeric(args[1])

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/dosage_functions.R'))

CHAIN = 'beta'

insertions = compile_mean_phenotype_data(c('v_gene', 'd_gene', 'd_gene', 'j_gene'), c('vd_insert', 'vd_insert', 'dj_insert', 'dj_insert'))

insertion_associations = compile_manhattan_plot_data(c('vd_insert', 'dj_insert'))

#validation cohort snp
snpID = 21722566 
genotypes = compile_all_genotypes_snp_list(snpID)

genotypes = association_genotype_assignment(-1, genotypes)

vd_boxplot = boxplot_by_snp(insertions, genotypes, 'vd_insert')
final_vd_insert_boxplot = vd_boxplot + xlab('rs3762093 SNP genotype') + ylab('Number of V-D-gene N-inserts')

dj_boxplot = boxplot_by_snp(insertions, genotypes, 'dj_insert')
final_dj_insert_boxplot = dj_boxplot + xlab('rs3762093 SNP genotype') + ylab('Number of D-J-gene N-inserts')

aligned_plots = align_plots(final_vd_insert_boxplot, final_dj_insert_boxplot, align = 'v')

together = plot_grid(aligned_plots[[1]], aligned_plots[[2]], nrow = 2, label_size = 40, labels = "AUTO", scale = 0.97)
ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/multipanel_insertion_boxplot_', CHAIN, '_rs3762093_discovery_cohort.pdf'), plot = together, width = 12, height = 17, units = 'in', dpi = 750, device = cairo_pdf)



