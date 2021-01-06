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
    } else if (phenotype_class == 'insertion_no_pca'){
        title = 'All N-insertion types \nwithout population structure correction'
    } else if (phenotype_class == 'p_addition_count' | phenotype_class == 'p_addition_fraction' | phenotype_class == 'p_addition_fraction_zero_trimming_subset'){
        title = 'All P-addition types conditioning on gene/cdr3'
    }
    return(title)
}

manhattan_plot <- function(dataframe, plot_title, file_name, plotting_p_value_cutoff, significance_cutoff = 0.05/35481497){
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
    
    snps = ggplot(dataframe) +
        geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = phenotype, shape = productive), alpha = 0.5, size = point_size) +
        geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.6) +
        facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") +
        theme_classic() +
        theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.text.x = element_blank()) +
        labs(y="-log10(p-value)", x="Chromosome Position") +
        geom_hline(yintercept=-log10(significance_cutoff), linetype="dashed", color = "green4", size=1.5) +
        ggtitle(plot_title) +
        scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8)) 
    
    ggsave(file_name, plot = snps, width = 25, height = 10, units = 'in', dpi = 750, device = 'pdf')
}
      
manhattan_plot_gene <- function(dataframe, plot_title, file_name, gene_subset){
    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

    gene_annotations = data.table(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)
    gene = gene_annotations[genes == gene_subset]
    dataframe = dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
    #Set plot parameters
    point_size = 5 
    significance_cutoff = 0.5/length(unique(dataframe$snp))
   
    snps = ggplot(dataframe) +
        geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = phenotype, shape = productive), alpha = 0.5, size = point_size) +
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
 
    
