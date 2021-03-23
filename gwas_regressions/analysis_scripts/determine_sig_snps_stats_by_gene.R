library(data.table)
setDTthreads(1)
library(plyr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))


trimming_data = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
trimming_naive_data = compile_manhattan_plot_data(c('v_trim_naive', 'j_trim_naive', 'd1_trim_naive', 'd0_trim_naive')) 
insertion_data = compile_manhattan_plot_data(c('vd_insert', 'dj_insert'))
pnuc_fraction_data = compile_manhattan_plot_data(c('v_pnucs_fraction_zero_trimming_subset', 'j_pnucs_fraction_zero_trimming_subset', 'd0_pnucs_fraction_zero_trimming_subset', 'd1_pnucs_fraction_zero_trimming_subset'))

sig_trim_artemis = get_sig_snps_stats(trimming_data, name = 'trimming', gene = 'artemis', 'genome-wide')
sig_trim_naive_artemis = get_sig_snps_stats(trimming_naive_data, name = 'trimming_naive', gene = 'artemis', 'genome-wide')
sig_trim_tcrb = get_sig_snps_stats(trimming_data, name = 'trimming', gene = 'tcrb', 'genome-wide')
sig_trim_naive_tcrb = get_sig_snps_stats(trimming_naive_data, name = 'trimming_naive', gene = 'tcrb', 'genome-wide')

sig_trim_naive_mhc = get_sig_snps_stats(trimming_naive_data, name = 'trimming_naive', gene = 'mhc', 'genome-wide')

sig_pnuc_frac_mhc = get_sig_snps_stats(pnuc_fraction_data, name = 'pnuc_fraction', gene = 'mhc', 'genome-wide')
sig_pnuc_frac_tcrb = get_sig_snps_stats(pnuc_fraction_data, name = 'pnuc_fraction', gene = 'tcrb', 'genome-wide')

sig_insert_dntt = get_sig_snps_stats(insertion_data, name = 'insertion', gene = 'dntt', 'genome-wide')
sig_insert_dntt_gene = get_sig_snps_stats(insertion_data, name = 'insertion', gene = 'dntt', 'gene-surround')


