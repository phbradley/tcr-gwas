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
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

ethnicity = fread(file = PCA_FILE)[,c('localID', 'scanID', 'race.g')]
pca = fread(PCA_FILE)
together = merge(ethnicity, pca, by = c('localID', 'scanID'))
setnames(together, 'race.g', 'ancestry_group')
together$pca_cluster = paste0('\"', together$ancestry_group, '\"-associated')

plot1 = ggparcoord(together, 
           columns = 4:35, 
           scale = 'uniminmax', 
           groupColumn = 'pca_cluster') +
geom_line(size=1.5, alpha = 0.4) +
theme_classic() + 
theme(text = element_text(size = 40)) +
    ggtitle('Population Structure PCA Scores by Ancestry Group') +
    xlab('\nPrincipal Component') +
    ylab('Scaled PC Scores\n') +
    scale_color_brewer(palette = 'Set2') 

final_plot1 = plot1 + theme_cowplot(font_family = 'Arial') + theme(legend.direction = "horizontal", text = element_text(size = 40), axis.text = element_text(size = 24), legend.text = element_text(size = 32), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(xlim = c(1,32), clip="off") + ggtitle('') + background_grid(major = 'xy') + labs(color = 'PCA Cluster') 

legend = get_legend(final_plot1 + theme(legend.position = c(0.39, 0.39)) + guides(color = guide_legend(nrow = 2)))

final_plot1_no_leg = final_plot1 + theme(legend.position = 'none')

scree_plot = data.frame(pc=1:32, varprop = pc$varprop)

plot2 = ggplot(scree_plot, aes(x = pc, y = 100*varprop)) +
    geom_line(size = 2) +
    geom_point(size = 6) +
    theme_classic() +
    theme(text = element_text(size = 40)) +
    ggtitle('Population Structure PCA Variance Explained') +
    xlab('\nPrincipal Component') +
    ylab('Proportion of Variance Explained\n') 

final_plot2 = plot2 + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 40), axis.text = element_text(size = 24), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(xlim = c(1,32), clip="off") + ggtitle('') + background_grid(major = 'xy')

aligned_plots = align_plots(final_plot1_no_leg, final_plot2, align = 'h', axis = 'bt')
together = plot_grid(aligned_plots[[2]], aligned_plots[[1]], ncol = 2, rel_widths = c(0.8, 2), label_size = 40, labels = "AUTO", scale = 0.97)

together_leg = plot_grid(together, legend, ncol = 1, rel_heights = c(2, 0.2), align = 'v', axis = 'l')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/pc_multipanel_plot.pdf'), plot = together_leg, width = 35, height = 13, units = 'in', dpi = 750, device = cairo_pdf)

