library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)
library(Cairo)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

dataframe = compile_manhattan_plot_data('gene_usage') 
significance_cutoff = 5.09e-11 

# get gene usage gwas data and associations by gene
tcrb = GENE_ANNOTATIONS[gene_common_name == 'tcrb']
mhc = GENE_ANNOTATIONS[gene_common_name == 'mhc']
znf = GENE_ANNOTATIONS[gene_common_name == 'znf443/znf709']

sig = dataframe[pvalue < significance_cutoff]
complement = fsetdiff(dataframe, sig) 
pairs = merge(sig[, -c('...1')], complement[, -c('...1')], by = c('snp', 'gene', 'parameter', 'phenotype', 'bootstraps', 'chr', 'hg19_pos'), all.x = TRUE)

# collapse double sig cases
collapse = sig
collapse[, count := .N, by = .(snp, gene)]
prod = collapse[count == 2 & productive == TRUE]
nonprod = collapse[count == 2 & productive == FALSE]
together = merge(prod[, -c('...1', 'count')], nonprod[, -c('...1', 'count')], by = c('snp', 'gene', 'parameter', 'phenotype', 'bootstraps', 'chr', 'hg19_pos'))

# get nonsig cases
sig_single = collapse[count == 1]
complement = fsetdiff(dataframe, sig_single[, -c('count')]) 
together_nonsig = merge(sig_single[, -c('...1', 'count')], complement[, -c('...1')], by = c('snp', 'gene', 'parameter', 'phenotype', 'bootstraps', 'chr', 'hg19_pos'), all.x = TRUE)

pairs = rbind(together, together_nonsig)

pairs[hg19_pos < (tcrb$pos2 + 200000) & hg19_pos > (tcrb$pos1 - 200000) & chr == tcrb$chr, region := 'TCRB']
pairs[hg19_pos < (mhc$pos2 + 200000) & hg19_pos > (mhc$pos1 - 200000) & chr == mhc$chr, region := 'MHC']
pairs[hg19_pos < (znf$pos2 + 200000) & hg19_pos > (znf$pos1 - 200000) & chr == znf$chr, region := 'ZNF443/ZNF709']
pairs[region == 'TCRB' & gene %in% TRB_NONPROD_ALLELE_GENES, region := 'TCRB (for usage of a gene with both\na productive and non-productive allele)']

pairs[productive.x == TRUE, prod_pvalue := pvalue.x]
pairs[productive.x == FALSE, nonprod_pvalue := pvalue.x]
pairs[productive.y == TRUE, prod_pvalue := pvalue.y]
pairs[productive.y == FALSE, nonprod_pvalue := pvalue.y]

pairs$region = factor(pairs$region, levels = c('MHC', 'TCRB', 'TCRB (for usage of a gene with both\na productive and non-productive allele)', 'ZNF443/ZNF709', 'Other'))
# randomize rows
set.seed(42)
rows = sample(nrow(pairs))
pairs = pairs[rows,]
pairs[is.na(region), region := 'Other']

plot = ggplot(pairs[region != 'Other']) +
    geom_point(aes(x = -log10(prod_pvalue), y = -log10(nonprod_pvalue), color = region), size = 4, alpha = 0.3) +
    geom_vline(xintercept = -log10(significance_cutoff), color = 'black', size = 2)+
    geom_hline(yintercept = -log10(significance_cutoff), color = 'black', size = 2)+
    geom_abline(intercept = 0, color = 'black', linetype = 'dashed', size = 2 )+
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 30, family = 'Arial'),axis.text.x=element_text(size = 18), axis.text.y = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    scale_color_discrete(name = "Locus") +
    ggtitle('') + 
    background_grid(major = 'xy') +
    xlab('productive -log10(p-value)') +
    ylab('non-productive -log10(p-value)') +
    panel_border(color = 'gray60', size = 1.5)+
    guides(color = guide_legend(override.aes = list(alpha=1)))

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/prod_nonprod_gene-usage_pvalue_comparison.pdf'), plot = plot, width = 18, height = 9, units = 'in', dpi = 750, device = cairo_pdf)


