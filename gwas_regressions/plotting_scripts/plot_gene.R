library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

# manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('tcr_div')), 
#                     phenotype_values = c('TCR diversity'),
#                     plot_title = set_manhattan_plot_title('tcr_div', 'mhc'), 
#                     file_name = make_file_name('tcr_div', 'mhc'), 
#                     gene_subset = 'mhc')

# manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('tcr_div')), 
#                     phenotype_values = c('TCR diversity'),
#                     plot_title = set_manhattan_plot_title('tcr_div', 'dntt'), 
#                     file_name = make_file_name('tcr_div', 'dntt'), 
#                     gene_subset = 'dntt')

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
                    phenotype_values = c('V-D-gene junction', 'D-J-gene junction'),
                    plot_title = set_manhattan_plot_title('insertion', 'dntt'), 
                    file_name = make_file_name('insertion', 'dntt'), 
                    gene_subset = 'dntt')

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('vd_insert_no_pca', 'dj_insert_no_pca')), 
                    phenotype_values = c('V-D-gene junction', 'D-J-gene junction'),
                    plot_title = set_manhattan_plot_title('insertion_no_pca', 'dntt'), 
                    file_name = make_file_name('insertion_no_pca', 'dntt'), 
                    gene_subset = 'dntt')

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
                    phenotype_values = c('V-gene trimming', 'J-gene trimming', '3\'-end D-gene trimming','5\'-end D-gene trimming'),
                    plot_title = set_manhattan_plot_title('trimming', 'tcrb'), 
                    file_name = make_file_name('trimming', 'tcrb'), 
                    gene_subset = 'tcrb', 
                    plotting_features = 'gene_features')

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
                    phenotype_values = c('V-gene trimming', 'J-gene trimming', '3\'-end D-gene trimming','5\'-end D-gene trimming'),
                    plot_title = set_manhattan_plot_title('trimming', 'tcrb'), 
                    file_name = make_file_name('trimming', 'tcrb_dj_zoom'), 
                    gene_subset = 'tcrb', 
                    plot_zoom = 'tcrb_dj_zoom', 
                    plotting_features = 'gene_features')

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
                    phenotype_values = c('V-gene trimming', 'J-gene trimming', '3\'-end D-gene trimming','5\'-end D-gene trimming'),
                    plot_title = set_manhattan_plot_title('trimming', 'tcrb'), 
                    file_name = make_file_name('trimming', 'tcrb_v_dj_zoom'), 
                    gene_subset = 'tcrb', 
                    plot_zoom = 'tcrb_v_dj_zoom', 
                    plotting_features = 'gene_features')




