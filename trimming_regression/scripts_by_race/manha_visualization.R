##library(plotly)
##library(manhattanly)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(GWASTools)
library(tidyverse)
#extrafont::loadfonts()
#omp_set_num_threads(1)
#blas_set_num_threads(1)
setDTthreads(threads = 1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/regression_functions.R')

# This script compiles all regression data from cluster across the entire genome
compile_by_race_data <- function(gene, trim_type, race_covariate, subset_by_race){
    stopifnot(race_covariate != subset_by_race)
    if (race_covariate == 'True'){
        filename = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/by_race/', gene, '_', trim_type, '_with_race_covariate.tsv')
    }
    if (subset_by_race == 'True'){
        filename = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/by_race/', gene, '_', trim_type, '_by_race.tsv')
    }
    data = fread(file = filename)[,-1]
    return(data)
}

find_significant_snps_by_gene <- function(trim_type, maf_data, random_effects, condensing, d_infer, repetitions,     maf_cutoff, pca_structure_correction, pca_type){
    maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)] 
    sig_data = compile_manhattan_plot_data(trim_type, maf_data, RANDOM_EFFECTS, CONDENSING, D_INFER, REPETITIONS, maf_cutoff, pca_structure_correction = 'False', pca_type = 'none')
}
compile_manhattan_plot_data <- function(trim_type, maf_data, random_effects, condensing, d_infer, repetitions, maf_cutoff, pca_structure_correction, pca_type){
    productive_file_name = make_compiled_regression_file_name(productivity = 'productive', trim_type,  pca_structure_correction, pca_type) 
     NOT_productive_file_name = make_compiled_regression_file_name(productivity = 'NOT_productive', trim_type, pca_structure_correction, pca_type) 
   
    productive_data = fread(productive_file_name, sep = "\t", fill=TRUE, header = TRUE)
    not_productive_data = fread(NOT_productive_file_name, sep = "\t", fill=TRUE, header = TRUE)
    
    data = rbind(productive_data, not_productive_data)[,-c(1,2)]
    data = merge(data, maf_data, by = 'snp')

    if (maf_cutoff != 'false'){
        data = data[maf > maf_cutoff]
    }

    data$trim_type = paste0(trim_type, '_', data$productivity)
    
    return(data)
}

# This script makes a manhattan plot for the entire genome!

manhattan_plot_cluster_by_race_gene <- function(trim_type, plotting_cutoff, subset, race_covariate, subset_by_race){
    bonferroni = 0.05/35481497
    
    data = compile_by_race_data(gene = subset, trim_type, race_covariate, subset_by_race)

    title = paste0(trim_type)
    if (race_covariate == 'True'){
        title = paste0(title, ' with race covariates')
        colnames(data) = c('snp', 'estimate', 'se', 'tvalue', 'pvalue', 'variable', 'chr', 'hg19_pos', 'productivity')
    }
    if (subset_by_race == 'True'){
        title = paste0(title, ' regressing by racial group')
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
    
    stopifnot(subset %in% c('False', 'artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra'))
    if (subset != 'False'){
        subset_gene = gene_annotations %>% filter(genes == subset)
        gene_annotations = subset_gene
        data = data %>% filter(chr == subset_gene$chr) %>% filter(hg19_pos < subset_gene$pos2 + 100000) %>% filter(hg19_pos > subset_gene$pos1 - 100000)
        point_size = 5
        alpha_gene = 0.1
    } else {
        point_size = 4 
        alpha_gene = 0.5
    }

     if (race_covariate == 'True'){
        data_subset = data[variable != 'snp' & variable != '(Intercept)']
        plot <- ggplot(data_subset) + 
            geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=variable), alpha = 0.5, size = point_size) + 
            geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = alpha_gene) + 
            facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + 
            theme_classic() + 
            #theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.title=element_text(size=25), legend.text=element_text(size=20), axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size= 20))+
            theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.text.x = element_text(size = 12, angle = 90))+ 
            labs(y="-log10(p-value)", x="Chromosome Position") + 
            ggtitle(title) + 
            scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8)) +
            guides(color = 'legend', shape = 'none', fill = 'legend')
    }
    if (subset_by_race == 'True'){
        data_subset = data[race != 'all_together' & race != 'Middle Eastern' & race != 'African']
        plot <- ggplot(data_subset) + 
            geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=race, shape=productivity), alpha = 0.5, size = point_size) + 
            geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = alpha_gene) + 
            facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + 
            theme_classic() + 
            #theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.title=element_text(size=25), legend.text=element_text(size=20), axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size= 20))+
            theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.text.x = element_text(size = 12, angle = 90))+ 
            labs(y="-log10(p-value)", x="Chromosome Position") + 
            ggtitle(title) + 
            scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8)) 
    }
    plot
}

