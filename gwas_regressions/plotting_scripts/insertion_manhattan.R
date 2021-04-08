library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

plotting_cutoff = -log10(5e-5)

features = c('vd_insert', 'dj_insert')
pretty_names = c('V-D-gene junction', 'D-J-gene junction')

manhattan_plot(dataframe = compile_manhattan_plot_data(features), 
               phenotype_values = pretty_names,     
               plot_title = set_manhattan_plot_title('insertion'), 
               file_name = make_file_name('insertion'), 
               plotting_p_value_cutoff = plotting_cutoff)
assign('insertion_manhattan', readRDS(paste0(make_file_name('insertion'), '.rds')))


revised_insertion = insertion_manhattan + theme_gray() + theme(panel.spacing.y=unit(3, "lines"), axis.ticks.x = element_blank(), panel.spacing.x=unit(0.25, "lines"), text = element_text(size = 40), axis.text.x = element_blank(), strip.text.x = element_text(size = 18)) + coord_cartesian(clip="off") + ggtitle('') + labs(color = 'Number of N-insertions')


file_name = paste0(make_file_name('insertion'), '.pdf')
ggsave(file_name, plot = revised_insertion, width = 25, height = 18, units= 'in', dpi = 750, device = 'pdf')

