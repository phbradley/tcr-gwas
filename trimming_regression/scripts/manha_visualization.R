##library(plotly)
##library(manhattanly)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(GWASTools)
library(tidyverse)
#omp_set_num_threads(1)
#blas_set_num_threads(1)
setDTthreads(threads = 1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/regression_functions.R')

# This script compiles all regression data from cluster across the entire genome

# This script compiles all minor allele fraction data (to be used in plotting cutoff)
compile_all_maf_data <- function(){
    data_files = list.files(path=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_results/'), pattern='*', full.names=TRUE)   

    print(paste0('compile file list'))
    assign(paste0('together'), NULL)
    count = 0
    for (file in data_files){
        # read file...
        if (file.size(file) == 1 | file.size(file) == 0){
            next
        }
        temp_file = fread(file, sep = "\t", fill=TRUE, header = TRUE)
        #together = rbind(together, temp_file)
        if (ncol(temp_file) > 2){
            assign(paste0('together'), rbindlist(list(get(paste0('together')), temp_file)))
        }
        count = count + 1
        print(paste0(count, ' of ', length(data_files), ' completed'))
    }
    
    file_name = paste0('maf_all_snps.tsv')

    write.table(together, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/', file_name), quote=FALSE, sep='\t', col.names = NA)
}

compile_manhattan_plot_data <- function(trim_type, maf_data, random_effects, condensing, bootstrap_count, maf_cutoff, bootstrap_rerun_count, pca_structure_correction){
    if (random_effects == 'True'){
         first_pass_file_name = paste0('_', trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_with_random_effects_', bootstrap_count, '_bootstraps.tsv')
         file_name = paste0('_', trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_with_random_effects_', bootstrap_rerun_count, '_bootstraps')
    } else {
         first_pass_file_name = paste0('_',trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_NO_random_effects_', bootstrap_count, '_bootstraps.tsv')
         file_name = paste0('_',trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_NO_random_effects_', bootstrap_rerun_count, '_bootstraps')
    }
    
    productive_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive', first_pass_file_name), sep = "\t", fill=TRUE, header = TRUE)
    not_productive_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive', first_pass_file_name), sep = "\t", fill=TRUE, header = TRUE)
    
    if (pca_structure_correction == 'True'){
        file_name = paste0(file_name, '_WITH_pca_structure_correction.tsv')
    } else{
        file_name = paste0(file_name, '.tsv')
    }

    prod_file_path = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive', file_name)
    NOTprod_file_path = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive', file_name)
        
    if (bootstrap_rerun_count != 'False' & file.exists(prod_file_path) & file.exists(NOTprod_file_path)){
        productive_boot_data = fread(prod_file_path, sep = "\t", fill=TRUE, header = TRUE)
        productive_boot_data = productive_boot_data[!duplicated(productive_boot_data$snp)]
        not_productive_boot_data = fread(NOTprod_file_path, sep = "\t", fill=TRUE, header = TRUE)
        not_productive_boot_data = not_productive_boot_data[!duplicated(not_productive_boot_data$snp)]
        data_boot = rbind(productive_boot_data, not_productive_boot_data)[,-c(1,2)]
        data_boot$bootstraps = paste0(bootstrap_rerun_count, '_bootstraps')

        # remove snps from non-bootstrapped dataframe if they have been boostrapped
        productive_data = productive_data[!productive_data$snp %in% productive_boot_data$snp,]
        not_productive_data = not_productive_data[!not_productive_data$snp %in% not_productive_boot_data$snp,]
        data = rbind(productive_data, not_productive_data)[,-c(1,2)]
        data$bootstraps = paste0(bootstrap_count, '_bootstraps')

        data = rbind(data, data_boot)
        data = merge(data, maf_data, by = 'snp')
    } else {
        data = rbind(productive_data, not_productive_data)[,-c(1,2)]
        data = merge(data, maf_data, by = 'snp')
        data$bootstraps = paste0(bootstrap_count, '_bootstraps')
    }

    if (maf_cutoff != 'False'){
        data = data[maf > maf_cutoff]
    }

    data$trim_type = paste0(trim_type, '_', data$productivity)
    
    return(data)
}

# This script makes a manhattan plot for the entire genome!

manhattan_plot_cluster <- function(trim_type, random_effects, bootstrap_count, condensing, plotting_cutoff, gene_annotations, maf_cutoff, bootstrap_rerun_count, pca_structure_correction){
    bonferroni = 0.05/35481497
    maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]

    if (trim_type == 'all_trim'){
        data = data.frame()
        for (trim in c('v_trim', 'd0_trim', 'd1_trim', 'j_trim')){
            data = rbind(data, compile_manhattan_plot_data(trim_type = trim, maf_data, random_effects, condensing, bootstrap_count, maf_cutoff, bootstrap_rerun_count, pca_structure_correction))
        }
    } else if (trim_type == 'all_insert'){
        data = data.frame()
        for (trim in c('vj_insert', 'vd_insert', 'dj_insert')){
            data = rbind(data, compile_manhattan_plot_data(trim_type = trim, maf_data, random_effects, condensing, bootstrap_count, maf_cutoff, bootstrap_rerun_count, pca_structure_correction))
        }
    } else {
        data = compile_manhattan_plot_data(trim_type, maf_data, random_effects, condensing, bootstrap_count, maf_cutoff, bootstrap_rerun_count, pca_structure_correction)
    }

    if (maf_cutoff != 'False'){
        title = paste0(trim_type, ' plot with ', bootstrap_count, ' bootstraps and ', maf_cutoff, ' MAF cutoff')
    } else {
        title = paste0(trim_type, ' plot with ', bootstrap_count, ' bootstraps')
    }

    if (pca_structure_correction == 'True'){
         title = paste0(title, ' with population structure PCA correction')
    }

    if (plotting_cutoff != 'False'){
        data_plotting_cutoff = data[-log10(pvalue) > plotting_cutoff]
        if (nrow(data)!=0){
            data = data_plotting_cutoff
        }
    }

    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

    gene_annotations = data.frame(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)

    if (gene_annotations == 'False'){
        ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=trim_type, shape = bootstraps), alpha = 0.5, size = 3) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
    } else {
        ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=trim_type, shape = bootstraps), alpha = 0.5, size = 3) + geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.5) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
    }
}

# This script makes a manhattan plot for a specific regression file

manhattan_plot_specific_file <- function(filename, trim_type, add_gene_annotations, maf_cutoff, bootstrap_count){
    bonferroni = 0.05/35481497
    data = fread(filename, sep = "\t", fill=TRUE, header = TRUE)
    maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]

    data = merge(data, maf_data, by = 'snp')

    if (maf_cutoff != 'False'){
        data = data[maf > maf_cutoff]
        title = paste0(trim_type, ' plot with ', bootstrap_count, ' bootstraps and ', maf_cutoff, ' MAF cutoff')
    } else {
        title = paste0(trim_type, ' plot with ', bootstrap_count, ' bootstraps')
    }

    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

    gene_annotations = data.frame(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)

    if (add_gene_annotations == 'False'){
        ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=productivity), alpha = 0.5, size = 2.5) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
    } else if (add_gene_annotations == 'True'){
        ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=productivity), alpha = 0.5, size = 2.5) + geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.5) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
    }
}
