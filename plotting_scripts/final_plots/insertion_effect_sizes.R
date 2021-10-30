library(cowplot)
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(Cairo)

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/analysis_functions.R'))

# get dntt gene level significant associations
sig_dntt = fread(paste0(OUTPUT_PATH, '/source_data/figure5-source-data2.txt'))

sig_dntt[productive == TRUE, productivity := 'productive']
sig_dntt[productive == FALSE, productivity := 'non-productive']

feature_colors = brewer.pal(n = max(4, length(unique(sig_dntt$phenotype))), name = "Set2")

plot = ggplot(sig_dntt, aes(y = effect_size, x = productivity)) +
    geom_boxplot(aes(fill = phenotype), lwd = 1.5) +
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
    facet_wrap(~phenotype)+
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 30, family = 'Arial'),legend.position = "none", axis.text.x=element_text(size = 18), axis.text.y = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    scale_fill_manual(guide = guide_legend(reverse=TRUE), values = feature_colors) +
    ggtitle('') + 
    background_grid(major = 'y') +
    ylab('Effect size (nucleotides)') +
    xlab('') +
    panel_border(color = 'gray60', size = 1.5)

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/gene_level_sig_dntt_effect_size_boxplot.pdf'), plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
