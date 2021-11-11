library(tidyverse)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(Cairo)
library(doParallel)
library(foreach)
library(rstatix)

args = commandArgs(trailingOnly=TRUE)

NCPU <<- as.numeric(args[1])

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/maf_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

# get insertion distribution and ethnicity data
insertions = compile_mean_phenotype_data(c('v_gene', 'd_gene', 'd_gene', 'j_gene'), c('vd_insert', 'vd_insert', 'dj_insert', 'dj_insert'))
ethnicity = fread(file = PCA_FILE)[,c('localID', 'race.g', 'race.s')]
together = merge(insertions, ethnicity, by = 'localID')[,c('localID', 'feature_of_interest', 'feature', 'race.g', 'productivity')]
average_all = together[, mean(feature_of_interest), by = .(productivity)]
setnames(average_all, 'V1', 'total_insert_mean')
setnames(together, 'race.g', 'ancestry_group')
together$ancestry_group = str_replace(together$ancestry_group, ' ', '\n')
together$pca_cluster = paste0('\"', together$ancestry_group, '\"-\nassociated')

# genotypes = compile_all_genotypes(as.numeric(10000), as.numeric(10))
# together = together[localID %in% genotypes$localID]

# t-test
t_test = together %>%
    group_by(productivity) %>%
    t_test(feature_of_interest ~ pca_cluster, ref.group = 'all')

outliers = together %>%
    group_by(productivity, pca_cluster) %>%
    identify_outliers(feature_of_interest) %>%
    as.data.table()

# print outliers from this analysis
print(outliers[, .N, by = .(localID, ancestry_group, race.s)])

plot = ggboxplot(together, x = 'pca_cluster', y = 'feature_of_interest', fill = 'pca_cluster', lwd = 1.5) +
    facet_grid(cols = vars(productivity)) +
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
    stat_pvalue_manual(t_test, label = "p = {p}", size = 8, y.position = 12, remove.bracket = TRUE, family = 'Arial') +
    theme_classic() + 
    theme(text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
    ggtitle('Mean Total Insertions by Ancestry Group') +
    geom_hline(data = average_all, aes(yintercept = total_insert_mean), size = 2.5, color = 'red', linetype = 2) +
    xlab('PCA cluster') +
    ylab('Mean total inserts') + 
    scale_fill_brewer(palette = 'Set2')

final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 40), axis.text.x=element_text(size = 22), axis.text.y = element_text(size = 22), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'xy') + panel_border(color = 'gray60', size = 1.5)
ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/insertions_by_race.pdf'), plot = final_plot, width = 25, height = 12, units = 'in', dpi = 750, device = cairo_pdf)

