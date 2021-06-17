library(data.table)
library(plyr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))


trimming_data = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
sig_trimming = get_association_count(trimming_data, name = 'trimming')

insertion_data = compile_manhattan_plot_data(c('vd_insert', 'dj_insert'))
sig_insertion = get_association_count(insertion_data, name = 'insertion')

pnuc_fraction_data = compile_manhattan_plot_data(c('v_pnucs_fraction_zero_trimming_subset', 'j_pnucs_fraction_zero_trimming_subset', 'd0_pnucs_fraction_zero_trimming_subset', 'd1_pnucs_fraction_zero_trimming_subset'))
sig_pnuc_fraction = get_association_count(pnuc_fraction_data, name = 'pnuc_fraction')

gene_usage_data = compile_manhattan_plot_data(c('gene_usage'))
sig_gene_usage = get_association_count(gene_usage_data, name = 'gene_usage')

