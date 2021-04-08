library(cowplot)
library(data.table)
setDTthreads(1)
library(plyr)
library(stringr)
library(RColorBrewer)
source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))

NCPU = 5

trimming_data = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = 'genome-wide', trimming_data)

tcrb = GENE_ANNOTATIONS[gene_common_name == 'tcrb']

tcrb_sig_snps = trimming_data[pvalue < significance_cutoff & chr == tcrb$chr & hg19_pos < tcrb$pos2 & hg19_pos > tcrb$pos1]

genotypes = compile_all_genotypes_snp_list(unique(tcrb_sig_snps$snp))
colnames(genotypes)[-1] = paste0('snp', colnames(genotypes)[-1])
d_gene_usage = get_d_gene_usage()

d_gene_2_alt_allele_usage = d_gene_usage[substring(d_gene, 1, 5) == 'TRBD2', alt_allele_prop := N/sum(N), by = localID][d_gene == 'TRBD2*02']


together = merge(genotypes, d_gene_2_alt_allele_usage[, c('localID', 'alt_allele_prop', 'd_gene')], by = 'localID')

for (snpID in colnames(genotypes)[-1]){
    boxplot_by_d_allele_usage(snpID, together, sig_snps = tcrb_sig_snps)
}

top_snp = paste0('snp', tcrb_sig_snps[order(pvalue)][1]$snp)

plot = boxplot_by_d_allele_usage(top_snp, together, sig_snps = tcrb_sig_snps, final_figure = TRUE)
filename = set_file_name(top_snp)

final_plot = plot +  theme_cowplot() + theme(legend.position = 'none', axis.text = element_text(size = 24), axis.line = element_blank(),text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) + ggtitle('') + background_grid(major = 'y') 

ggsave(filename, plot = final_plot, width = 14, height = 10, units = 'in', dpi = 750, device = 'pdf')

compare = list(c('0', '1'), c('1', '2'), c('0', '2'))
compare_means( as.formula(paste0('alt_allele_prop ~', top_snp)), comparisons = compare, p.adjust.method = "bonferroni", method='t.test', data = together)
