library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('tcr_div')), 
               phenotype_values = c('TCR diversity'),
               plot_title = set_manhattan_plot_title('tcr_div'), 
               file_name = make_file_name('tcr_div'), 
               plotting_p_value_cutoff = -log10(5e-4))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_pnucs_count', 'j_pnucs_count', 'd0_pnucs_count', 'd1_pnucs_count')), 
               phenotype_values = c('V-gene P-nucleotide count', 'J-gene P-nucleotide count', '5\'-end D-gene P-nucleotide count', '3\'-end D-gene P-nucleotide count'), 
               plot_title = set_manhattan_plot_title('p_addition_count'), 
               file_name = make_file_name('p_addition_count'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_pnucs_fraction', 'j_pnucs_fraction', 'd0_pnucs_fraction', 'd1_pnucs_fraction')), 
               phenotype_values = c('V-gene P-nucleotide TCR fraction', 'J-gene P-nucleotide TCR fraction', '5\'-end D-gene P-nucleotide TCR fraction', '3\'-end D-gene P-nucleotide TCR fraction'),
               plot_title = set_manhattan_plot_title('p_addition_fraction'), 
               file_name = make_file_name('p_addition_fraction'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim')), 
               phenotype_values = c('V-gene trimming', 'J-gene trimming', '3\'-end D-gene trimming', '5\'-end D-gene trimming'), 
               plot_title = set_manhattan_plot_title('trimming'), 
               file_name = make_file_name('trimming'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_trim_naive', 'j_trim_naive', 'd1_trim_naive', 'd0_trim_naive')), 
               phenotype_values = c('V-gene trimming', 'J-gene trimming', '3\'-end D-gene trimming','5\'-end D-gene trimming'),
               plot_title = set_manhattan_plot_title('trimming_naive'), 
               file_name = make_file_name('trimming_naive'), 
               plotting_p_value_cutoff = -log10(5e-5))


manhattan_plot(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
               phenotype_values = c('V-D N-inserts', 'D-J N-inserts'),
               plot_title = set_manhattan_plot_title('insertion'), 
               file_name = make_file_name('insertion'), 
               plotting_p_value_cutoff = -log10(5e-5), 
               significance_cutoff_type = c('genome-wide', 'dntt'))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('vd_insert_no_pca', 'dj_insert_no_pca')), 
               phenotype_values = c('V-D N-inserts', 'D-J N-inserts'),
               plot_title = set_manhattan_plot_title('insertion_no_pca'), 
               file_name = make_file_name('insertion_no_pca'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_pnucs_fraction_zero_trimming_subset', 'j_pnucs_fraction_zero_trimming_subset', 'd0_pnucs_fraction_zero_trimming_subset', 'd1_pnucs_fraction_zero_trimming_subset')), 
               phenotype_values = c('V-gene P-nucleotide non-trimmed TCR fraction', 'J-gene P-nucleotide non-trimmed TCR fraction', '5\'-end D-gene P-nucleotide non-trimmed TCR fraction', '3\'-end D-gene P-nucleotide non-trimmed TCR fraction'),
               plot_title = set_manhattan_plot_title('p_addition_fraction_zero_trimming_subset'), 
               file_name = make_file_name('p_addition_fraction_zero_trimming_subset'), 
               plotting_p_value_cutoff = -log10(5e-5))

manhattan_plot(dataframe = compile_manhattan_plot_data(c('v_trim_zero_trimming_fraction', 'j_trim_zero_trimming_fraction', 'd1_trim_zero_trimming_fraction', 'd0_trim_zero_trimming_fraction')), 
               phenotype_values = c('non-V-gene-trimmed TCR fraction', 'non-J-gene-trimmed TCR fraction', 'non-3\'-D-gene-trimmed TCR fraction', 'non-5\'-D-gene-trimmed TCR fraction'),
               plot_title = set_manhattan_plot_title('zero_trimming_fraction'), 
               file_name = make_file_name('zero_trimming_fraction'), 
               plotting_p_value_cutoff = -log10(5e-5))
