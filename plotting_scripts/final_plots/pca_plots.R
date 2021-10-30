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
library(stringr)
library(tidyverse)

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

together = fread(PCA_FILE)
setnames(together, 'race.g', 'ancestry_group')
together$pca_cluster = paste0('\"', together$ancestry_group, '\"-associated')
colnames(together) = str_remove_all(colnames(together), "EV")
#scatter rows: 
set.seed(42)
rows = sample(nrow(together))
together = together[rows,]

plot1 = ggparcoord(together, 
           columns = 3:34, 
           scale = 'uniminmax', 
           groupColumn = 'pca_cluster') +
geom_line(size=3, alpha = 0.6) +
theme_classic() + 
theme(text = element_text(size = 40)) +
    ggtitle('Population Structure PCA Scores by Ancestry Group') +
    xlab('\nPrincipal Component') +
    ylab('Scaled PC Scores\n') +
    scale_color_brewer(palette = 'Set2') 

final_plot1 = plot1 + theme_cowplot(font_family = 'Arial') + theme(legend.direction = "horizontal", text = element_text(size = 40), axis.text = element_text(size = 24), legend.text = element_text(size = 32), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(xlim = c(1,32), clip="off") + ggtitle('') + background_grid(major = 'xy') + labs(color = 'PCA Cluster') 

legend = get_legend(final_plot1 + theme(legend.position = c(0.25, 0.39)) + guides(color = guide_legend(nrow = 2, override.aes = list(size = 7))))

final_plot1_no_leg = final_plot1 + theme(legend.position = 'none')

scree_plot = fread(PCA_VARIANCE_FILE)

scaleFUN <- function(x) sprintf("%.2f", x)

plot2 = ggplot(scree_plot, aes(x = as.factor(pc), y = 100*variance_explained, group=1)) +
    geom_line(size = 2, stat='summary', fun.y=sum) +
    stat_summary(fun.y=sum, geom="line")+
    geom_point(size = 6) +
    scale_y_continuous(labels=scaleFUN) +
    theme_classic() +
    theme(text = element_text(size = 40)) +
    ggtitle('Population Structure PCA Variance Explained') +
    xlab('\nPrincipal Component') +
    ylab('Proportion of\nVariance Explained\n')

final_plot2 = plot2 + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 40), axis.text = element_text(size = 24), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + ggtitle('') + background_grid(major = 'xy')

aligned_plots = align_plots(final_plot1_no_leg, final_plot2, align = 'v', axis = 'lr')
grid_together = plot_grid(aligned_plots[[2]], aligned_plots[[1]], nrow = 2, rel_heights = c(1.5, 4), label_size = 40, labels = "AUTO", scale = 0.97)

together_leg = plot_grid(grid_together, legend, ncol = 1, rel_heights = c(2, 0.2), align = 'v', axis = 'l')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/pc_multipanel_plot.pdf'), plot = together_leg, width = 35, height = 25, units = 'in', dpi = 750, device = cairo_pdf)

