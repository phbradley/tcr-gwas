library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)
library(Cairo)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

# get gene usage gwas data and associations by gene
gene = GENE_ANNOTATIONS[gene_common_name == 'mhc']

pos1 = gene$pos1
pos2 = gene$pos2
chrom = gene$chr

dataframe = fread(paste0(OUTPUT_PATH, '/source_data/figure1-source-data1.txt'))
# filter for mhc region associations
mhc_associations = dataframe[hg19_pos < (pos2 + 200000) & hg19_pos > (pos1 - 200000) & chr == chrom]

significance_cutoff = 5.09e-11 
sig_mhc = mhc_associations[pvalue < significance_cutoff]

# make sure effect sizes are going in the correct direction
med = sig_mhc[, median(effect_size), by = .(productive, gene)]
setnames(med, 'V1', 'med_effect_size')

q1 = sig_mhc[, quantile(effect_size, .25), by = .(productive, gene)]
setnames(q1, 'V1', 'q1_effect_size')

q3 = sig_mhc[, quantile(effect_size, .75), by = .(productive, gene)]
setnames(q3, 'V1', 'q3_effect_size')
together = merge(med, q1)
together = merge(together, q3)

together$gene = factor(together$gene, levels = unique(together[order(med_effect_size)]$gene))

plot = ggplot() +
    geom_vline(xintercept = 0, color = 'gray60', size = 2)+
    geom_point(data = together, aes(x = med_effect_size, y = gene), size = 8) +
    geom_segment(data = together, aes(x = q1_effect_size, xend = q3_effect_size, y = gene, yend = gene), size = 3)+
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 30, family = 'Arial'),legend.position = "none", axis.text.x=element_text(size = 18), axis.text.y = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    ggtitle('') + 
    background_grid(major = 'y') +
    xlab('Effect size (proportion of repertoire using gene)') +
    ylab('') +
    panel_border(color = 'gray60', size = 1.5)

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/sig_mhc_effect_size_boxplot.pdf'), plot = plot, width = 14, height = 8, units = 'in', dpi = 750, device = cairo_pdf)
#
