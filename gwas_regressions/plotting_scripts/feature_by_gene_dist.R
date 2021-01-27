library(tidyverse)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
#TODO functions to create these data
insertions = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/insertions_by_patient.tsv')[,-1]
setnames(insertions, 'patient_id', 'localID')
trimmings = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/trim_by_patient.tsv')[,-1]
setnames(trimmings, 'patient_id', 'localID')

create_distribution_data <- function(data, gene_types, paired_feature_types){
    stopifnot(length(gene_types)==length(paired_feature_types))
    together = data.table()
    for (i in seq(length(gene_types))){
        temp = data
        names(temp)[names(temp) == paste0(gene_types[i])] <- 'gene_of_interest'
        names(temp)[names(temp) == paste0(paired_feature_types[i])] <- 'feature_of_interest'
        temp$productivity = ifelse(data$productive == FALSE, 'non-productive', 'productive')
        temp_condensed = temp[, mean(feature_of_interest), by = .(gene_of_interest, productivity, localID)]
        temp_condensed$feature = paired_feature_types[i]
        temp_condensed$gene_type = gene_types[i]
        together = rbind(together, temp_condensed)
        temp = data.table()
        temp_condensed = data.table()
    }
    setnames(together, 'V1', 'feature_of_interest')
    return(together)
}

create_title <- function(paired_feature_types){
    type = strsplit(paired_feature_types, '_')[[1]][2]
    if (type == 'trim'){
        title = 'Trimming distribution by gene'
    } else if (type == 'insert'){
        title = 'N-insertion distribution by gene'
    } else if (type == 'pnucs'){
        title = 'P-addition fraction distribution by gene'
    }
    return(title)
}

create_filename <- function(paired_feature_types){
    type = strsplit(paired_feature_types, '_')[[1]][2]
    if (type == 'trim'){
        title = 'trimming_dist_by_gene.pdf'
    } else if (type == 'insert'){
        title = 'n-insertion_dist_by_gene.pdf'
    } else if (type == 'pnucs'){
        title = 'p-addition_dist_by_gene.pdf'
    }
    return(title)
}


plot_distributions_by_gene <- function(data, gene_types, paired_feature_types){
    data_condensed = create_distribution_data(data, gene_types, paired_feature_types)
    data_condensed = data_condensed[gene_of_interest != '-']
    data_condensed$facet_variable = paste0(data_condensed$feature, ' by ', data_condensed$gene_type)
    mean_dt = data_condensed %>%
        group_by(gene_of_interest, productivity, facet_variable) %>%
        summarize(median= median(feature_of_interest))
    plot = ggplot(data_condensed, aes(x=feature_of_interest)) +
        facet_grid(cols = vars(productivity), rows = vars(facet_variable)) +
        geom_vline(data = mean_dt, aes(xintercept=median, color = gene_of_interest), alpha = 0.1, linetype='solid', size=1) +
        geom_density(aes(color = gene_of_interest), alpha = 0.3, adjust = 2.5) + 
        theme_classic() +
        theme(text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
        ggtitle(create_title(paired_feature_types))
    filename = create_filename(paired_feature_types)
    ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/', filename), plot = plot, width = 15, height = 20, units = 'in', dpi = 750, device = 'pdf')
}

plot_distributions_by_gene(insertions, c('v_gene', 'd_gene', 'd_gene', 'j_gene'), c('vd_insert', 'vd_insert', 'dj_insert', 'dj_insert'))
plot_distributions_by_gene(trimmings, c('v_gene', 'd_gene', 'd_gene', 'j_gene'), c('v_trim', 'd0_trim', 'd1_trim', 'j_trim'))

#TODO FINISH this!!
