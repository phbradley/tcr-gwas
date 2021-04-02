library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

plotting_cutoff = -log10(5e-5)

CORRECTION_TYPE <<- 'allele_status' 
stopifnot(CORRECTION_TYPE %in% c('snp_interaction', 'snp_allele_status_interaction', 'allele_status'))

# PHENOTYPE_TYPE <<- 'p-addition_fraction_trimming_subset'
PHENOTYPE_TYPE <<- 'trimming'
stopifnot(PHENOTYPE_TYPE %in% c('trimming', 'p-addition_fraction_trimming_subset'))

interaction_corrected_data = fread(paste0(OUTPUT_PATH, '/results/d_allele_linkage/', PHENOTYPE_TYPE, '_regressions_', CORRECTION_TYPE, '_correction_by_gene_cdr3_d_infer-True_8_PCAir_PCs_tcrb.tsv')) 

dataframe_correction = interaction_corrected_data[linkage_correction == 'With correction']
dataframe_no_correction = interaction_corrected_data[linkage_correction == 'Without correction']

manhattan_plot_gene_compare(dataframe_correction,  
                            dataframe_no_correction, 
                            # phenotype_values = c('V-gene', '5\'-end D-gene', '3\'-end D-gene', 'J-gene'),
                            phenotype_values = c('5\'-end D-gene'),
                            plot_title = NULL, 
                            file_name = get_file_name_compare('tcrb'),
                            gene_subset = 'tcrb')


