library(tidyverse)
library(GWASTools)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(GGally)
library(Cairo)

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

ethnicity = fread(file = PCA_FILE)[,c('localID', 'scanID', 'race.g')]

counts = ethnicity[, .N, by = race.g]
counts$race.g = factor(counts$race.g, levels = c('Caucasian', 'Hispanic', 'African', 'Asian', 'Native American', 'Middle Eastern'))
counts$race.g = str_replace(counts$race.g, ' ', '\n')
counts$pca_cluster = paste0('\"', counts$race.g, '\"-\nassociated')

plot = ggplot(counts) +
    geom_bar(aes(x = pca_cluster, y = N, fill = pca_cluster), stat = 'identity') +
    xlab('PCA cluster') +
    ylab('Subject count') +
    scale_fill_brewer(palette = 'Set2')

final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = 'none', text = element_text(size = 22), axis.text = element_text(size = 12), axis.text.x = element_text(size = 12), axis.line = element_blank(), axis.ticks = element_blank()) + background_grid(major = 'y') 

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/population_composition.pdf'), plot = final_plot, width = 8, height = 8, units = 'in', dpi = 750, device = cairo_pdf)


