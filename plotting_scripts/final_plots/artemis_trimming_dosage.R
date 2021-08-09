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

# compile trimming distribution data
trimmings = compile_mean_phenotype_data(c('v_gene', 'd_gene', 'd_gene', 'j_gene'), c('v_trim', 'd0_trim', 'd1_trim', 'j_trim'))

# compile trimming gwas data
trimming_associations = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
gene = GENE_ANNOTATIONS[gene_common_name == 'artemis']

# filter for artemis region associations
artemis_associations = trimming_associations[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
top_associations = artemis_associations[order(pvalue)][1:10]
top_associations[, min_p := min(pvalue), by = .(phenotype, productive)]
top_j_snp = top_associations[phenotype == 'j_trim'][1]
top_v_snp = top_associations[phenotype == 'v_trim'][1]

# get genotypes for top snps
j_genotypes = compile_all_genotypes_snp_list(top_j_snp$snp)
v_genotypes = compile_all_genotypes_snp_list(top_v_snp$snp)
j_genotypes = association_genotype_assignment(top_j_snp$slope, j_genotypes)
v_genotypes = association_genotype_assignment(top_v_snp$slope, v_genotypes)

# get gene usage
v_gene_usage = trimmings[, .N, by = .(v_gene, localID, productive)]
v_gene_usage[, max_N := max(N), by = .(localID, productive)]
top_v_gene = v_gene_usage[max_N == N][, .N, by = .(v_gene)][order(-N)][1]
j_gene_usage = trimmings[, .N, by = .(j_gene, localID, productive)]
j_gene_usage[, max_N := max(N), by = .(localID, productive)]

#use the top TRBJ1 gene (so that there is no potential for D-gene misidentification biases)
top_j_gene = j_gene_usage[max_N == N][, .N, by = .(j_gene)][order(-N)][2]


j_trim_boxplot = boxplot_by_snp(trimmings, j_genotypes,'j_trim', top_j_gene$j_gene)
final_j_trim_boxplot = j_trim_boxplot + xlab('rs41298872 SNP genotype') + ylab('Number of J-gene nucleotides deleted')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/j_trim_boxplot_', top_j_gene$j_gene, '_snp', top_j_snp$snp, '.pdf'), plot = final_j_trim_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

v_trim_boxplot = boxplot_by_snp(trimmings, v_genotypes,'v_trim', top_v_gene$v_gene)
final_v_trim_boxplot = v_trim_boxplot + xlab('Genotype of the top V-gene trimming associated DCLRE1C SNP') + ylab('Number of V-gene nucleotides deleted')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/v_trim_boxplot_', top_v_gene$v_gene, '_snp', top_v_snp$snp, '.pdf'), plot = final_v_trim_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)


