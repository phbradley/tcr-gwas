library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)

source('config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plot_src/plotting_functions/manhattan_plot_functions.R'))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_pnucs_count', 'j_pnucs_count', 'd0_pnucs_count', 'd1_pnucs_count')), 
               plot_title = set_manhattan_plot_title('p_addition_count'), 
               file_name = make_file_name('p_addition_count'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_pnucs_fraction', 'j_pnucs_fraction', 'd0_pnucs_fraction', 'd1_pnucs_fraction')), 
               plot_title = set_manhattan_plot_title('p_addition_fraction'), 
               file_name = make_file_name('p_addition_fraction'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
               plot_title = set_manhattan_plot_title('trimming'), 
               file_name = make_file_name('trimming'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_trim_naive', 'j_trim_naive', 'd1_trim_naive', 'd0_trim_naive')), 
               plot_title = set_manhattan_plot_title('trimming_naive'), 
               file_name = make_file_name('trimming_naive'), 
               plotting_p_value_cutoff = -log10(5e-5))


manhattan_plot(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
               plot_title = set_manhattan_plot_title('insertion'), 
               file_name = make_file_name('insertion'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('total_insert')), 
               plot_title = set_manhattan_plot_title('insertion'), 
               file_name = make_file_name('insertion'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('vd_insert_no_pca', 'dj_insert_no_pca')), 
               plot_title = set_manhattan_plot_title('insertion_no_pca'), 
               file_name = make_file_name('insertion_no_pca'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_pnucs_fraction_zero_trimming_subset', 'j_pnucs_fraction_zero_trimming_subset', 'd0_pnucs_fraction_zero_trimming_subset', 'd1_pnucs_fraction_zero_trimming_subset')), 
               plot_title = set_manhattan_plot_title('p_addition_fraction_zero_trimming_subset'), 
               file_name = make_file_name('p_addition_fraction_zero_trimming_subset'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_trim_zero_trimming_fraction', 'j_trim_zero_trimming_fraction', 'd1_trim_zero_trimming_fraction', 'd0_trim_zero_trimming_fraction')), 
               plot_title = set_manhattan_plot_title('zero_trimming_fraction'), 
               file_name = make_file_name('zero_trimming_fraction'), 
               plotting_p_value_cutoff = -log10(5e-5))
