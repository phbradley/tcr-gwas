library(GWASTools)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(GGally)
library(Cairo)

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))[,c('localID', 'scanID', 'race.g')]
pc = getobj(paste0(PROJECT_PATH, "/tcr-gwas/_ignore/pcair_pcair.RData"))
pc_scores = as.data.frame(pc$vectors)
colnames(pc_scores) = as.character(seq(1:32))
pc_scores$scanID = as.integer(rownames(pc_scores))
together = merge(ethnicity, pc_scores, by = 'scanID')
setnames(together, 'race.g', 'ancestry_group')

plot1 = ggparcoord(together, 
           columns = 4:35, 
           scale = 'uniminmax', 
           groupColumn = 'ancestry_group') +
geom_line(size=1.5, alpha = 0.4) +
theme_classic() + 
theme(text = element_text(size = 40)) +
    ggtitle('Population Structure PCA Scores by Ancestry Group') +
    xlab('Principal Component') +
    ylab('Scaled PC Scores') +
    scale_color_brewer(palette = 'Set2')

final_plot1 = plot1 + theme_cowplot(font_family = 'Arial') + theme(text = element_text(size = 40), axis.text = element_text(size = 24), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'xy') + labs(color = 'Ancestry group')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/pc_parallel_coordinates_by_race.pdf'), plot = final_plot1, width = 30, height = 12, units = 'in', dpi = 750, device = cairo_pdf)

scree_plot = data.frame(pc=1:32, varprop = pc$varprop)

plot2 = ggplot(scree_plot, aes(x = pc, y = 100*varprop)) +
    geom_line(size = 2) +
    geom_point(size = 6) +
    theme_classic() +
    theme(text = element_text(size = 40)) +
    ggtitle('Population Structure PCA Variance Explained') +
    xlab('Principal Component') +
    ylab('Proportion of Variance Explained')

final_plot2 = plot2 + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 40), axis.text = element_text(size = 24), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'xy') 

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/pc_scree_plot.pdf'), plot = final_plot2, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
