library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

plotting_cutoff = -log10(5e-5)

features = c('gene_usage')

manhattan_plot(dataframe = compile_manhattan_plot_data(features), 
               plot_title = set_manhattan_plot_title('gene_usage'), 
               file_name = make_file_name('gene_usage'), 
               plotting_p_value_cutoff = plotting_cutoff, 
               gene_usage = TRUE)

assign('gene_usage_manhattan', readRDS(paste0(make_file_name('gene_usage'), '.rds')))

revised_gene_usage = gene_usage_manhattan + theme_cowplot() + theme(panel.spacing.y=unit(3, "lines"), axis.ticks.x = element_blank(), panel.spacing.x=unit(0, "lines"), strip.placement = 'outside', strip.background.x = element_blank(), text = element_text(size = 40), axis.text.x = element_blank(), strip.text.x = element_text(size = 12), axis.text.y = element_text(size = 24), axis.line.y = element_line(color = 'gray60', size = 1.5), axis.line.x = element_blank(), axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + labs(color = 'Gene Usage') + geom_hline(aes(yintercept = Inf), size = 1.5, color = 'gray60')+geom_hline(aes(yintercept = -Inf), size = 1.5, color = 'gray60') + scale_x_continuous() + geom_point(aes(x = Inf, y = -Inf), size = 3, color = 'gray60', shape = '\u007c')

file_name = paste0(make_file_name('gene_usage'), '.pdf')
ggsave(file_name, plot = revised_gene_usage + background_grid(major = 'y'), width = 25, height = 18, units= 'in', dpi = 750, device = 'pdf')


