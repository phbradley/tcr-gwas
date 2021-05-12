library(data.table)
library(plyr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))


trimming_data = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
insertion_data = compile_manhattan_plot_data(c('vd_insert', 'dj_insert'))

find_feature_overlaps(trimming_data, gene_subset = 'artemis')
find_feature_overlaps(trimming_data, gene_subset = 'tcrb')
find_feature_overlaps(insertion_data, gene_subset = 'dntt', significance_cutoff_type = 'dntt')

