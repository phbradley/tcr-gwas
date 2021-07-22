library(data.table)
library(plyr)

source('config/config.R')
source('config/file_paths.R')

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))


trimming_data = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
trimming_lambdas = get_lambdas(trimming_data, 'trimming')

trimming_naive_data = compile_manhattan_plot_data(c('v_trim_naive', 'j_trim_naive', 'd1_trim_naive', 'd0_trim_naive')) 
trimming_naive_lambdas = get_lambdas(trimming_naive_data, 'trimming_naive')

insertion_data = compile_manhattan_plot_data(c('vd_insert', 'dj_insert'))
insertion_lambdas = get_lambdas(insertion_data, 'insertion')

pnuc_fraction_data = compile_manhattan_plot_data(c('v_pnucs_fraction_zero_trimming_subset', 'j_pnucs_fraction_zero_trimming_subset', 'd0_pnucs_fraction_zero_trimming_subset', 'd1_pnucs_fraction_zero_trimming_subset'))
pnuc_fraction_lambdas = get_lambdas(pnuc_fraction_data, 'pnucs_fraction_zero_trimming_subset')

gene_usage_data = compile_manhattan_plot_data(c('gene_usage'))
gene_usage_lambdas = get_lambdas_gene_usage(gene_usage_data)
