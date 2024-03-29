library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)
library(ggh4x)
library(ggnewscale)


source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

manhattan_plot_loci(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
                    phenotype_values = c('V-gene\ntrimming', 'J-gene\ntrimming', '3\'-end D-gene\ntrimming', '5\'-end D-gene\ntrimming'),
                    plot_title = set_manhattan_plot_title('trimming', 'artemis'), 
                    file_name = make_file_name('trimming', 'artemis_loci'), 
                    gene_subset = 'artemis')

assign('artemis_trim_loci', readRDS(paste0(make_file_name('trimming', 'artemis_loci'), '.rds')))

file_name = paste0(make_file_name('trimming', 'artemis_loci'), '.pdf')
final_plot = artemis_trim_loci + ggtitle('') + theme_cowplot(font_family = 'Arial') + theme(text = element_text(size = 40), axis.text.y = element_text(size = 30), axis.line = element_blank(), axis.text.x = element_blank()) + coord_cartesian(clip="off") + background_grid(major = 'xy') + xlab('DCLRE1C position') + ylab('') 

ggsave(file_name, plot= final_plot, width = 35, height = 10, units= 'in', dpi = 750, device = cairo_pdf)
