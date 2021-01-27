source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/regression_parameters.R'))

get_compiled_file_name_with_arguments <- function(phenotype, condensing_variable, infer_missing_d_gene, pca_count){
    file_name = paste0(OUTPUT_PATH, '/results/', phenotype, '_regressions_', condensing_variable, '_d_infer-', infer_missing_d_gene, '_', pca_count, '_PCAir_PCs.tsv')
    return(file_name)
}

make_file_name <- function(phenotype_class, gene_subset = 'NA'){
    file = paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/all_', phenotype_class, '_genome_wide_manhattan.pdf')
    if (gene_subset != 'NA'){
        file = paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/all_', phenotype_class, '_', gene_subset, '_manhattan.pdf')
    }
    return(file)
}


compile_manhattan_plot_data <- function(phenotype_list){
    compiled_data = data.table()
    for (phenotype in phenotype_list){
        parameters = set_regression_parameters(phenotype)
        file_name = get_compiled_file_name_with_arguments(parameters[['phenotype']], parameters[['condensing_variable']], parameters[['infer_missing_d_gene']], parameters[['pca_count']])
        compiled_data = rbind(compiled_data, fread(file_name))
    }
    return(compiled_data)
}

set_manhattan_plot_title <- function(phenotype_class, gene_subset = 'NA'){
    if (phenotype_class == 'trimming' | phenotype_class == 'zero_trimming_fraction'){
        title = 'All trimming types conditioning on gene/cdr3'
    } else if (phenotype_class == 'insertion'){
        title = 'All N-insertion types'
    } else if (phenotype_class == 'total_insertion'){
        title = 'Total N-insertion'
    } else if (phenotype_class == 'trimming_naive'){
        title = 'All trimming types without gene/cdr3 conditioning'
        if (gene_subset != 'NA'){
            title = 'All trimming types\nwithout gene/cdr3 conditioning'
        }
    } else if (phenotype_class == 'insertion_no_pca'){
        title = 'All N-insertion types \nwithout population structure correction'
    } else if (phenotype_class == 'p_addition_count' | phenotype_class == 'p_addition_fraction' | phenotype_class == 'p_addition_fraction_zero_trimming_subset'){
        title = 'All P-addition types conditioning on gene/cdr3'
        if (gene_subset != 'NA'){
            title = 'All P-addition types\nconditioning on gene/cdr3'
        }
    } else if (phenotype_class == 'tcr_div'){
        title = 'TCR Diversity'
    }
    return(title)
}

determine_significance_cutoff <- function(alpha, type, genome_wide_dataframe){
    stopifnot(type %in% c('genome-wide', 'artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra'))
    if (type == 'genome-wide'){
        sig_cutoff = alpha/35481497
    } else {
        genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
        chr = c(10, 6, 10, 11, 7, 14)
        pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
        pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

        gene_annotations = data.table(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)
        gene = gene_annotations[genes == type]
        dataframe = genome_wide_dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
        sig_cutoff = alpha/length(unique(dataframe$snp))
    }
    return(sig_cutoff)
}

manhattan_plot <- function(dataframe, phenotype_values = NA, plot_title, file_name, plotting_p_value_cutoff, significance_cutoff_type = 'genome-wide'){
    stopifnot(length(significance_cutoff_type) <= 2)
    if (!is.na(phenotype_values)){
        dataframe$feature = mapvalues(dataframe$phenotype, from = unique(dataframe$phenotype), to = phenotype_values) 
    } else {
        dataframe$feature = dataframe$phenotype
    }
    
    significance_cutoff = c()
    sig_colors = c()
    if (length(significance_cutoff_type) == 2){
        significance_cutoff = c(significance_cutoff, determine_significance_cutoff(0.05, significance_cutoff_type[1], dataframe))
        significance_cutoff = c(significance_cutoff, determine_significance_cutoff(0.05, significance_cutoff_type[2], dataframe))
        sig_colors = c(sig_colors, 'green4')
        sig_colors = c(sig_colors, 'dodgerblue2')
    } else {
        significance_cutoff = c(significance_cutoff, determine_significance_cutoff(0.05, significance_cutoff_type, dataframe))
        sig_colors = c(sig_colors, 'green4')
    }
    sig_lines = data.frame(significance = paste0(significance_cutoff_type, ' significance cutoff'), cutoff = significance_cutoff)
    names(sig_colors) = sig_lines$significance

    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)
    
    gene_annotations = data.frame(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)
    gene_annotations$pos1 = gene_annotations$pos1 - 1000000
    gene_annotations$pos2 = gene_annotations$pos2 + 1000000
    dataframe = dataframe[-log10(pvalue)>plotting_p_value_cutoff]
    #Set plot parameters
    point_size = 5 
    
    if (!('productive' %in% colnames(dataframe))) {
        dataframe$productivity = 'all'
    } else {
        dataframe[productive == TRUE, productivity := 'productive']
        dataframe[productive == FALSE, productivity := 'non-productive']
    }
    
    feature_colors = brewer.pal(n = max(3, length(unique(dataframe$feature))), name = "Set2")
    if (length(feature_colors) > length(unique(dataframe$feature))){
        feature_colors = feature_colors[1:length(unique(dataframe$feature))]
    }
    names(feature_colors) = unique(dataframe$feature)

    snps = ggplot(dataframe) +
        geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = feature, shape = productivity), alpha = 0.5, size = point_size) +
        geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.6) +
        facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") +
        theme_classic() +
        theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.text.x = element_blank()) +
        labs(y="-log10(p-value)", x="Chromosome Position") +
        geom_hline(data = sig_lines[1,], aes(yintercept=-log10(cutoff), linetype=significance), color = 'green4', size=2) +
        ggtitle(plot_title) +
        scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
    if (nrow(sig_lines) > 1){
        snps = snps + 
            geom_hline(data = sig_lines[2,], aes(yintercept=-log10(cutoff), linetype=significance), color = 'dodgerblue2', size=2) +
            scale_linetype_manual(name = 'significance', values = c(2,2), guide = guide_legend(override.aes = list(color = c('green4', 'dodgerblue2'))))
    } else {
        snps = snps +
            scale_linetype_manual(name = 'significance', values = 2, guide = guide_legend(override.aes = list(color = 'green4')))
    }
    ggsave(file_name, plot = snps, width = 25, height = 12, units = 'in', dpi = 750, device = 'pdf')
}
      
manhattan_plot_gene <- function(dataframe, phenotype_values = NA, plot_title, file_name, gene_subset){
    if (!is.na(phenotype_values)){
        dataframe$feature = mapvalues(dataframe$phenotype, from = unique(dataframe$phenotype), to = phenotype_values) 
    } else {
        dataframe$feature = dataframe$phenotype
    }
    significance_cutoff = determine_significance_cutoff(0.05, significance_cutoff_type = gene_subset, dataframe)
    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

    gene_annotations = data.table(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)
    gene = gene_annotations[genes == gene_subset]
    dataframe = dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
    #Set plot parameters
    point_size = 5 
   
    if (!('productive' %in% colnames(dataframe))) {
        dataframe$productivity = 'all'
    } else {
        dataframe[productive == TRUE, productivity := 'productive']
        dataframe[productive == FALSE, productivity := 'non-productive']
    }

    snps = ggplot(dataframe) +
        geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = feature, shape = productivity), alpha = 0.5, size = point_size) +
        geom_rect(data = gene, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.1) +
        facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") +
        theme_classic() +
        theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.text.x = element_blank()) +
        labs(y="-log10(p-value)", x="Chromosome Position") +
        geom_hline(yintercept=-log10(significance_cutoff), linetype="dashed", color = "green4", size=1.5) +
        ggtitle(plot_title) +
        scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8)) 
    
    ggsave(file_name, plot = snps, width = 13, height = 10, units = 'in', dpi = 750, device = 'pdf')
}
 
manhattan_plot_loci <- function(dataframe, phenotype_values = NA, plot_title, file_name, gene_subset, significance_cutoff_type = 'genome-wide'){
    if (!is.na(phenotype_values)){
        dataframe$feature = mapvalues(dataframe$phenotype, from = unique(dataframe$phenotype), to = phenotype_values) 
    } else {
        dataframe$feature = dataframe$phenotype
    }

    significance_cutoff = determine_significance_cutoff(0.05, significance_cutoff_type, dataframe)

    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

    gene_annotations = data.table(gene_locus = genes, chr = chr, pos1 = pos1, pos2 = pos2)
    gene = gene_annotations[gene_locus == gene_subset]
    if (significance_cutoff == 'zoom'){
        dataframe = dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
        significance_cutoff = 0.5/length(unique(dataframe$snp))
        plot_title = paste0(plot_title, '\nGene level Bonferroni significance cutoff')
    } else {
        dataframe = dataframe[hg19_pos < (gene$pos2) & hg19_pos > (gene$pos1) & chr == gene$chr]
        plot_title = paste0(plot_title, '\nGenome-wide Bonferroni significance cutoff')
    }
    
    sig_snps = dataframe[pvalue < significance_cutoff]

    if (!('productive' %in% colnames(dataframe))) {
        dataframe$productive = 'all'
    }
    sig_snps$productivity = ifelse(sig_snps$productive == FALSE, 'non-productive', 'productive')
    sig_snps$parameter = paste0(sig_snps$feature, '_', sig_snps$productivity_name)
    most_sig_snps = sig_snps[, .SD[which.min(pvalue)], by = .(feature, productivity)] 
    most_sig_snps$most_significant_association = "most significant association"
    snps = ggplot(sig_snps) +
        geom_point(aes(x = hg19_pos, y = productivity, size = -log10(pvalue)), alpha = 0.25) +
        geom_point(data = most_sig_snps, aes(x = hg19_pos, y = productivity, color = most_significant_association, size = -log10(pvalue)), color = 'dodgerblue', stroke = 4, shape = 1) +
        scale_size_continuous(range = c(4, 40)) +
        geom_rect(data = gene, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = gene_locus), alpha = 0.08) +
        facet_wrap(~feature, ncol = 1, strip.position = 'left') +
        theme_classic() +
        theme(text = element_text(size = 40), axis.text.x = element_blank()) +
        labs(x="Position") +
        ggtitle(plot_title)
    
    ggsave(file_name, plot = snps, width = 35, height = 10, units = 'in', dpi = 750, device = 'pdf')
}
    
