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

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/dosage_functions.R'))

CHAIN <<- 'beta'

trimmings = fread(file = MEAN_INSERTS)[,-1]
setnames(trimmings, 'patient_id', 'localID')

trimming_associations = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
gene = GENE_ANNOTATIONS[gene_common_name == 'artemis']

# validation cohort snp for trimming
snpID = 20717849

genotypes = compile_all_genotypes_snp_list(snpID)
genotypes = association_genotype_assignment(trimming_associations[snp == snpID]$slope, genotypes)

v_gene_usage = trimmings[, .N, by = .(v_gene, localID, productive)]
v_gene_usage[, max_N := max(N), by = .(localID, productive)]

top_v_gene = v_gene_usage[max_N == N][, .N, by = .(v_gene)][order(-N)][1]

j_gene_usage = trimmings[, .N, by = .(j_gene, localID, productive)]
j_gene_usage[, max_N := max(N), by = .(localID, productive)]

#use the top TRBJ1 gene (so that there is no potential for D-gene misidentification biases)
top_j_gene = j_gene_usage[max_N == N][, .N, by = .(j_gene)][order(-N)][2]


j_trim_boxplot = boxplot_by_snp(trimmings, genotypes, 'j_trim', top_j_gene$j_gene)
final_j_trim_boxplot = j_trim_boxplot + xlab('rs12768894 SNP genotype') + ylab('Number of J-gene nucleotides deleted')

v_trim_boxplot = boxplot_by_snp(trimmings, genotypes, 'v_trim', top_v_gene$v_gene)
final_v_trim_boxplot = v_trim_boxplot + xlab('rs12768894 SNP genotype') + ylab('Number of V-gene nucleotides deleted')

aligned_plots = align_plots(final_v_trim_boxplot, final_j_trim_boxplot, align = 'v')

together = plot_grid(aligned_plots[[1]], aligned_plots[[2]], nrow = 2, label_size = 40, labels = "AUTO", scale = 0.97)
ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/multipanel_trimming_boxplot_', CHAIN, '_rs12768894_discovery_cohort.pdf'), plot = together, width = 12, height = 17, units = 'in', dpi = 750, device = cairo_pdf)



