library(cowplot)
library(data.table)
setDTthreads(1)
library(plyr)
library(stringr)
library(RColorBrewer)
library(Cairo)
library(heatmaply)

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))

NCPU = 5

# compile trimming gwas data
trimming_data = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))

# get significance cutoff
significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = 'genome-wide', trimming_data)

# get tcrb associations
tcrb = GENE_ANNOTATIONS[gene_common_name == 'tcrb']
tcrb_sig_snps = trimming_data[pvalue < significance_cutoff & chr == tcrb$chr & hg19_pos < tcrb$pos2 & hg19_pos > tcrb$pos1]

# get genotypes
genotypes = compile_all_genotypes_snp_list(unique(tcrb_sig_snps$snp))
colnames(genotypes)[-1] = paste0('snp', colnames(genotypes)[-1])

# get d gene alleles
allele_statuses = fread(D_ALLELES)
allele_statuses[allele_0 == 'TRBD2*02', alt_allele_genotype := 1]
allele_statuses[allele_0 != 'TRBD2*02', alt_allele_genotype := 0]
allele_statuses[allele_1 == 'TRBD2*02', alt_allele_genotype := alt_allele_genotype + 1]
colnames(allele_statuses) = c('allele_0', 'allele_1', 'localID', 'TRBD2_alt_allele_genotype')

# merge data and get top snp
together = merge(genotypes, allele_statuses, by = 'localID')
top_snp = paste0('snp', tcrb_sig_snps[order(pvalue)][1]$snp)

# D allele boxplot
plot = boxplot_by_d_allele_genotype(top_snp, together, sig_snps = tcrb_sig_snps, final_figure = TRUE)
cols = c('TRBD2_alt_allele_genotype', top_snp)
subset = together[,..cols]

# chi squared test
chi = chisq.test(table(subset))
print(chi)

# reformat data
heatmap_data = subset[!is.na(snp16814596), .N, by = .(TRBD2_alt_allele_genotype , snp16814596)]

plot = ggplot(heatmap_data, aes(x = TRBD2_alt_allele_genotype, snp16814596))+ 
    geom_tile(aes(fill = N)) +
    scale_fill_gradient(low = 'white', high = "royalblue4", name = 'Subject count') 

filename = set_file_name(top_snp, type = 'allele_genotype', final_figure = TRUE)

final_plot = plot +  theme_cowplot(font_family = 'Arial') + theme(axis.line = element_blank(), axis.text = element_text(size = 20), text = element_text(size = 30), legend.key.width = unit(0.5, "in"), legend.key.height = unit(1, "in")) + ggtitle('') + ylab('rs2367486 genotype') + xlab('TRBD2*02 allele genotype') +
panel_border(color = 'gray60', size = 1.5)


ggsave(filename, plot = final_plot, width = 12, height = 10, units = 'in', dpi = 750, device = cairo_pdf)


