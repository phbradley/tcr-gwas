library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
               plot_title = set_manhattan_plot_title('insertion', 'dntt'), 
               file_name = make_file_name('insertion', 'dntt'), 
               gene_subset = 'dntt')

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('vd_insert_no_pca', 'dj_insert_no_pca')), 
               plot_title = set_manhattan_plot_title('insertion_no_pca', 'dntt'), 
               file_name = make_file_name('insertion_no_pca', 'dntt'), 
               gene_subset = 'dntt')

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('v_trim', 'd0_trim', 'd1_trim', 'j_trim')), 
               plot_title = set_manhattan_plot_title('trimming', 'artemis'), 
               file_name = make_file_name('trimming', 'artemis'), 
               gene_subset = 'artemis')

