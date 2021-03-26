library(data.table)
setDTthreads(1)
library(plyr)
library(xtable)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))


trimming_data = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
trimming_lambdas = get_lambdas_xtable(trimming_data, 'trimming')

trimming_naive_data = compile_manhattan_plot_data(c('v_trim_naive', 'j_trim_naive', 'd1_trim_naive', 'd0_trim_naive')) 
trimming_naive_lambdas = get_lambdas_xtable(trimming_naive_data, 'trimming_naive')

insertion_data = compile_manhattan_plot_data(c('vd_insert', 'dj_insert'))
insertion_lambdas = get_lambdas_xtable(insertion_data, 'insertion')

pnuc_fraction_data = compile_manhattan_plot_data(c('v_pnucs_fraction_zero_trimming_subset', 'j_pnucs_fraction_zero_trimming_subset', 'd0_pnucs_fraction_zero_trimming_subset', 'd1_pnucs_fraction_zero_trimming_subset'))
pnuc_fraction_lambdas = get_lambdas_xtable(pnuc_fraction_data, 'pnucs_fraction_zero_trimming_subset')

# gene_usage_data = compile_manhattan_plot_data(c('gene_usage'))
# gene_usage_lambdas = get_lambdas_xtable_gene_usage(gene_usage_data)
# sig_gene_usage_tcrb = get_sig_snps_stats(gene_usage_data, name = 'gene_usage', gene = 'tcrb', 'genome-wide')
# sig_gene_usage_mhc = get_sig_snps_stats(gene_usage_data, name = 'gene_usage', gene = 'mhc', 'genome-wide')

sig_trim_artemis = get_sig_snps_stats(trimming_data, name = 'trimming', gene = 'artemis', 'genome-wide')
sig_trim_naive_artemis = get_sig_snps_stats(trimming_naive_data, name = 'trimming_naive', gene = 'artemis', 'genome-wide')
sig_trim_tcrb = get_sig_snps_stats(trimming_data, name = 'trimming', gene = 'tcrb', 'genome-wide')
sig_trim_naive_tcrb = get_sig_snps_stats(trimming_naive_data, name = 'trimming_naive', gene = 'tcrb', 'genome-wide')

sig_trim_naive_mhc = get_sig_snps_stats(trimming_naive_data, name = 'trimming_naive', gene = 'mhc', 'genome-wide')

sig_pnuc_frac_mhc = get_sig_snps_stats(pnuc_fraction_data, name = 'pnuc_fraction', gene = 'mhc', 'genome-wide')
sig_pnuc_frac_tcrb = get_sig_snps_stats(pnuc_fraction_data, name = 'pnuc_fraction', gene = 'tcrb', 'genome-wide')

sig_insert_dntt = get_sig_snps_stats(insertion_data, name = 'insertion', gene = 'dntt', 'genome-wide')
sig_insert_dntt_gene = get_sig_snps_stats(insertion_data, name = 'insertion', gene = 'dntt', 'gene-surround')

# Explore gene usage at tcrb more...
genes = unique(sig_gene_usage_tcrb$gene)
prod_v_gene_percent nrow(sig_gene_usage_tcrb[substring(gene, 4, 4) == 'V'][, .N, by = .(gene, productive)][productive== TRUE])/length(genes[substring(genes, 4, 4) == 'V'])
not_prod_v_gene_percent nrow(sig_gene_usage_tcrb[substring(gene, 4, 4) == 'V'][, .N, by = .(gene, productive)][productive== FALSE])/length(genes[substring(genes, 4, 4) == 'V'])


#chr19 sig peak

chr19_peak = gene_usage_data[chr == 19 & pvalue < 4.72e-11]
