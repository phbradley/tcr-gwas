library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)
library(ggh4x)
library(ggnewscale)
library(Cairo)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

manhattan_plot_loci(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
                    phenotype_values = c('V-D-gene\nN-inserts', 'D-J-gene\nN-inserts'),
                    plot_title = set_manhattan_plot_title('insertion', 'dntt'), 
                    file_name = make_file_name('insertion', 'dntt_loci_zoom'), 
                    gene_subset = 'dntt',
                    significance_cutoff = 'dntt')

assign('dntt_insert_loci', readRDS(paste0(make_file_name('insertion', 'dntt_loci_zoom'), '.rds')))

file_name = paste0(make_file_name('insertion', 'dntt_loci_zoom'), '.pdf')

final_plot = dntt_insert_loci + ggtitle('') + theme_cowplot(font_family = 'Arial') + theme(text = element_text(size = 40), axis.text.y = element_text(size = 30),axis.line = element_blank(), axis.text.x = element_blank()) + coord_cartesian(clip="off") + background_grid(major = 'xy') + xlab('Extended DNTT locus position') + ylab('') 

ggsave(file_name, plot= final_plot, width = 35, height = 10, units= 'in', dpi = 750, device = cairo_pdf)
