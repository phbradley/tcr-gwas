library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)

source('config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plot_src/plotting_functions/manhattan_plot_functions.R'))

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
               plot_title = set_manhattan_plot_title('insertion', 'dntt'), 
               file_name = make_file_name('insertion', 'dntt'), 
               gene_subset = 'dntt')

