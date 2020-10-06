library(data.table)
library(ggplot2)
library(RColorBrewer)
library(GWASTools)
library(tidyverse)


manhattan_plot_cluster <- function(trim_type, random_effects, bootstrap_count, plotting_cutoff, gene_annotations, maf_cutoff, bootstrap_rerun_count){
    if (random_effects == 'True'){
         file_name = paste0('_', trim_type, '_snps_regression_with_weighting_condensing_by_gene_with_random_effects_', bootstrap_count, '_bootstraps.tsv')
    } else {
         file_name = paste0('_',trim_type, '_snps_regression_with_weighting_condensing_by_gene_NO_random_effects_', bootstrap_count, '_bootstraps.tsv')
    }
    productive_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive', file_name), sep = "\t", fill=TRUE, header = TRUE)[,-c(1,2)]
    not_productive_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive', file_name), sep = "\t", fill=TRUE, header = TRUE)[,-c(1,2)]

    if (random_effects == 'True'){
        file_name = paste0('_', trim_type, '_snps_regression_with_weighting_condensing_by_gene_with_random_effects_', bootstrap_rerun_count, '_bootstraps.tsv')
    } else {
        file_name = paste0('_',trim_type, '_snps_regression_with_weighting_condensing_by_gene_NO_random_effects_', bootstrap_rerun_count, '_bootstraps.tsv')
    }

    productive_boot_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive', file_name), sep = "\t", fill=TRUE, header = TRUE)[,-c(1,2)]
    colnames(productive_boot_data) = c('snp', 'chr', 'hg_19_pos', 'intercept_boostrapped', 'slope_boostrapped', 'standard_error_boostrapped', 'pvalue_bootstrapped', 'productivity')
    
    not_productive_boot_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive', file_name), sep = "\t", fill=TRUE, header = TRUE)[,-c(1,2)]
    colnames(not_productive_boot_data) = c('snp', 'chr', 'hg_19_pos', 'intercept_boostrapped', 'slope_boostrapped', 'standard_error_boostrapped', 'pvalue_bootstrapped', 'productivity')
    

    # remove snps from non-bootstrapped dataframe if they have been boostrapped
    productive_data_pre = productive_data[productive_data$snp %in% productive_boot_data$snp,]
    not_productive_data_pre = not_productive_data[not_productive_data$snp %in% not_productive_boot_data$snp,]

    together_productive = merge(productive_data_pre, productive_boot_data)
    together_productive = together_productive[!duplicated(together_productive$snp)]

    together_not_productive = merge(not_productive_data_pre, not_productive_boot_data)
    together_not_productive = together_not_productive[!duplicated(together_not_productive$snp)]

    together = rbind(together_productive,  together_not_productive)


