library(cowplot)
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(Cairo)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
#TODO functions to create these data
insertions = fread(file = MEAN_INSERTS)[,-1]
setnames(insertions, 'patient_id', 'localID')
trimmings = fread(file = MEAN_TRIMS)[,-1]
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


plot_distributions_by_gene <- function(data, gene_types, paired_feature_types, plot_variable_names, pretty_gene_names, xlab){
    data_condensed = create_distribution_data(data, gene_types, paired_feature_types)
    data_condensed = data_condensed[gene_of_interest != '-']
    data_condensed$facet_variable = mapvalues(data_condensed$feature, from = unique(data_condensed$feature), to = plot_variable_names)
    data_condensed$pretty_gene = mapvalues(data_condensed$gene_type, from = unique(data_condensed$gene_type), to = pretty_gene_names)
    data_condensed$facet_variable = paste0(data_condensed$facet_variable, '\nby ', data_condensed$pretty_gene)

    mean_dt = data_condensed[, median(feature_of_interest), by = .(gene_of_interest, productivity, facet_variable)]
    colnames(mean_dt) = c('gene_of_interest', 'productivity', 'facet_variable', 'median')
    
    plot = ggplot(data_condensed, aes(x=feature_of_interest, color = facet_variable, group = gene_of_interest)) +
        facet_grid(cols = vars(productivity), rows = vars(facet_variable)) +
        # geom_vline(data = mean_dt, aes(xintercept=median, color = facet_variable), alpha = 0.1, linetype='solid', size=1) +
        # geom_density(aes(group = gene_of_interest, color = facet_variable), alpha = 0.3, adjust = 2.5) + 
        stat_ecdf(geom = "step", size = 1, alpha = 0.5) +
        # geom_rug(sides="b", color = 'gray40', alpha = 0.5) +
        theme_classic() +
        theme(text = element_text(size = 40), legend.position = "none") +
        ggtitle(create_title(paired_feature_types))+
        scale_color_brewer(palette = 'Set2')+
        xlim(0, 25)
    
    final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", axis.text = element_text(size = 24), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'xy') + panel_border(color = 'gray60', size = 1.5) + xlab(xlab) + ylab('Cumulative probability')

    filename = create_filename(paired_feature_types)
    ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/', filename), plot = final_plot, width = 15, height = 20, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_distributions_by_gene(insertions, c('v_gene', 'd_gene', 'd_gene', 'j_gene'), c('vd_insert', 'vd_insert', 'dj_insert', 'dj_insert'), c('V-D-gene insertions', 'D-J-gene insertions'), c('V-gene','D-gene', 'J-gene'), xlab = 'Number of N-insertions')

plot_distributions_by_gene(trimmings, c('v_gene', 'd_gene', 'd_gene', 'j_gene'), c('v_trim', 'd0_trim', 'd1_trim', 'j_trim'), c('V-gene trimming', '5\'-D-gene trimming', '3\'-end-D-gene trimming', 'J-gene trimming'), c('V-gene', 'D-gene', 'J-gene'), xlab = 'Number of trimmed nucleotides')


#TODO FINISH this!!