library(data.table)
setDTthreads(1)
library(plyr)
library(stringr)
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
