library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggh4x)
library(ggnewscale)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

manhattan_plot_loci(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
                    phenotype_values = c('V-gene\ntrimming', 'J-gene\ntrimming', '3\'-end D-gene\ntrimming', '5\'-end D-gene\ntrimming'),
                    plot_title = set_manhattan_plot_title('trimming', 'artemis'), 
                    file_name = make_file_name('trimming', 'artemis_loci'), 
                    gene_subset = 'artemis')

manhattan_plot_loci(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
                    phenotype_values = c('V-gene\ntrimming', 'J-gene\ntrimming', '3\'-end D-gene\ntrimming','5\'-end D-gene\ntrimming'),
                    plot_title = set_manhattan_plot_title('trimming', 'tcrb'), 
                    file_name = make_file_name('trimming', 'tcrb_loci'), 
                    gene_subset = 'tcrb')

manhattan_plot_loci(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
                    phenotype_values = c('V-gene\ntrimming', 'J-gene\ntrimming', '3\'-end D-gene\ntrimming','5\'-end D-gene\ntrimming'),
                    plot_title = set_manhattan_plot_title('trimming', 'tcrb'), 
                    file_name = make_file_name('trimming', 'tcrb_loci_zoom'), 
                    gene_subset = 'tcrb', 
                    significance_cutoff = 'tcrb', 
                    plot_zoom = 'tcrb_dj_zoom')



manhattan_plot_loci(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
                    phenotype_values = c('V-D-gene\nN-inserts', 'D-J-gene\nN-inserts'),
                    plot_title = set_manhattan_plot_title('insertion', 'dntt'), 
                    file_name = make_file_name('insertion', 'dntt_loci'), 
                    gene_subset = 'dntt')

manhattan_plot_loci(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
                    phenotype_values = c('V-D-gene\nN-inserts', 'D-J-gene\nN-inserts'),
                    plot_title = set_manhattan_plot_title('insertion', 'dntt'), 
                    file_name = make_file_name('insertion', 'dntt_loci_zoom'), 
                    gene_subset = 'dntt',
                    significance_cutoff = 'dntt')

