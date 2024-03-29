source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/regression_parameters.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/gene_annotations.R'))

get_compiled_file_name_with_arguments <- function(phenotype, condensing_variable, keep_missing_d_gene, pca_count){
    file_name = paste0(OUTPUT_PATH, '/results/', phenotype, '_regressions_', condensing_variable, '_d_infer-', keep_missing_d_gene, '_', pca_count, '_PCAir_PCs.tsv')
    return(file_name)
}

make_file_name <- function(phenotype_class, gene_subset = 'NA', cohort = 'NA'){
    file = paste0(PROJECT_PATH, '/tcr-gwas/figures/all_', phenotype_class, '_genome_wide_manhattan')
    if (gene_subset != 'NA'){
        file = paste0(PROJECT_PATH, '/tcr-gwas/figures/all_', phenotype_class, '_', gene_subset, '_manhattan')
    }
    if (cohort != 'NA') {
        file = paste0(file, '_', cohort)
    }
    return(file)
}

compile_manhattan_plot_data <- function(phenotype_list){
    compiled_data = data.table()
    for (phenotype in phenotype_list){
        parameters = set_regression_parameters(phenotype)
        file_name = get_compiled_file_name_with_arguments(parameters[['phenotype']], parameters[['condensing_variable']], parameters[['keep_missing_d_gene']], parameters[['pca_count']])
        compiled_data = rbind(compiled_data, fread(file_name))
    }
    qc_file = fread(SNP_META_DATA_FILE)
    not_pass = qc_file[quality.filter == FALSE | quality.filter.cc == FALSE]
    filtered_snps = setdiff(compiled_data$snp, not_pass$snp)
    compiled_data = compiled_data[snp %in% filtered_snps]
    return(compiled_data)
}

extract_subject_ID <- function(tcr_repertoire_file_path){
    # This function extracts the subjectID from the file name
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][8]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

create_distribution_data <- function(data, gene_types, paired_feature_types){
    stopifnot(length(gene_types)==length(paired_feature_types))
    type = str_split(paired_feature_types[1], '_')[[1]][2]
    together = data.table()
    print(type)
    if (type == 'trim'){
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
    } else if (type == 'insert'){
        together = data
        together$productivity = ifelse(data$productive == FALSE, 'non-productive', 'productive')
        inserts = unique(paired_feature_types)
        together$feature_of_interest = 0
        for (insert in inserts){
            together$feature_of_interest = together$feature_of_interest + together[[insert]]
        }
        together = together[, mean(feature_of_interest), by = .(productivity, localID)]
        together$feature = 'total_insert'
    }
    setnames(together, 'V1', 'feature_of_interest')
    return(together)
}

compile_mean_phenotype_data <- function(gene_types, paired_feature_types){
    files = list.files(TCR_REPERTOIRE_DATA_DIRECTORY, pattern = "*.tsv", full.names=TRUE)
    registerDoParallel(cores=NCPU)
    mean_tcr_data = foreach(file = files, .combine = 'rbind') %dopar% {
        file_data = fread(file)
        file_data$localID = extract_subject_ID(file)
        create_distribution_data(file_data, gene_types, paired_feature_types)
    }
    stopImplicitCluster()
    return(mean_tcr_data)
}

set_manhattan_plot_title <- function(phenotype_class, gene_subset = 'NA'){
    if (phenotype_class == 'trimming'){
        title = 'Amount of gene-specific nucleotide trimming,\nconditioning on gene/cdr3'
    }
    else if (phenotype_class == 'zero_trimming_fraction'){
        title = 'Fraction of non-gene-trimmed TCRs, conditioning on gene/cdr3'
    } else if (phenotype_class == 'insertion'){
        title = 'Number of N-insertions'
    } else if (phenotype_class == 'total_insertion'){
        title = 'Total number of N-insertions'
    } else if (phenotype_class == 'trimming_naive'){
        title = 'Amount of gene-specific nucleotide trimming without gene/cdr3 conditioning'
        if (gene_subset != 'NA'){
            title = 'Amount of gene-specific nucleotide trimming\nwithout gene/cdr3 conditioning'
        }
    } else if (phenotype_class == 'insertion_no_pca'){
        title = 'Number of N-insertions\nwithout population structure correction'
    } else if (phenotype_class == 'p_addition_count'){
        title = 'Number of gene-specific P-nucleotides, conditioning on gene/cdr3'
    } else if (phenotype_class == 'p_addition_fraction'){
        title = 'Fraction of TCRs containing P-nucleotides, conditioning on gene/cdr3'
    } else if (phenotype_class == 'p_addition_fraction_zero_trimming_subset'){
        title = 'Fraction of un-trimmed TCRs containing P-nucleotides, \nconditioning on gene/cdr3'
        # if (gene_subset != 'NA'){
        #     title = 'All P-addition types\nconditioning on gene/cdr3'
        # }
    } else if (phenotype_class == 'tcr_div'){
        title = 'TCR Diversity'
    } else if (phenotype_class == 'gene_usage'){
        title = 'Gene usage'
    } else if (phenotype_class == 'gene_usage_no_PCA'){
        title = 'Gene usage without population structure correction'
    } else if (phenotype_class == 'missing_d_fraction'){
        title = 'Fraction of TCRs containing a \"missing D gene\"'
    } else {
        title = ''
    } 
    return(title)
}

determine_significance_cutoff <- function(alpha, type, genome_wide_dataframe, phenotype = unique(genome_wide_dataframe$phenotype)){
    stopifnot(type %in% c('genome-wide', 'artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra'))
    features = length(phenotype)*length(unique(genome_wide_dataframe$productive))
    if (type == 'genome-wide'){
        snps_total = length(unique(genome_wide_dataframe$snp))
    } else {
        gene = GENE_ANNOTATIONS[gene_common_name == type]
        gene_chr = gene$chr
        gene_pos1 = gene$pos1
        gene_pos2 = gene$pos2
 
        dataframe = genome_wide_dataframe[hg19_pos < (gene_pos2 + 200000) & hg19_pos > (gene_pos1 - 200000) & chr == gene_chr]
        snps_total = length(unique(dataframe$snp))
    }
    sig_cutoff = (alpha/snps_total)/features
    return(sig_cutoff)
}

manhattan_plot <- function(dataframe, phenotype_subset = NA, phenotype_values = NA, plot_title, file_name, plotting_p_value_cutoff, significance_cutoff_type = 'genome-wide', gene_usage = NA, subsample_point_cutoff = NA){
    stopifnot(length(significance_cutoff_type) <= 2)
    dataframe = dataframe[!is.na(phenotype)]
    if (!is.na(phenotype_values)){
        dataframe$feature = mapvalues(dataframe$phenotype, from = unique(dataframe$phenotype), to = phenotype_values) 
    } else {
        if (!is.na(gene_usage)){
            dataframe$feature = paste0(substring(dataframe$gene, 4, 4), '-gene usage')
        } else{
            dataframe$feature = dataframe$phenotype
        }
    }
    
    if (!is.na(gene_usage)){
        significance_cutoff = determine_significance_cutoff(0.05, significance_cutoff_type, dataframe, phenotype = unique(dataframe$gene))
    } else {
        significance_cutoff = determine_significance_cutoff(0.05, significance_cutoff_type, dataframe)
    }

    sig_lines = data.frame(significance = paste0(significance_cutoff_type, ' significance cutoff'), cutoff = significance_cutoff)
 
    # set phenotype colors
    dataframe = dataframe[order(-feature)]
    feature_colors = brewer.pal(n = max(4, length(unique(dataframe$feature))), name = "Set2")
    names(feature_colors) = unique(dataframe$feature)

   
    GENE_ANNOTATIONS$mid_pos = rowMeans(GENE_ANNOTATIONS[, c('pos1', 'pos2')])
    dataframe = dataframe[-log10(pvalue)>plotting_p_value_cutoff]
    if (!is.na(phenotype_subset)){
        dataframe = dataframe[phenotype %in% phenotype_subset]
    }

    feature_colors = feature_colors[names(feature_colors) %in% unique(dataframe$feature)]
    #Set plot parameters
    point_size = 5 
    
    if (!('productive' %in% colnames(dataframe))) {
        dataframe$productivity = 'all'
    } else {
        dataframe[productive == TRUE, productivity := 'productive']
        dataframe[productive == FALSE, productivity := 'non-productive']
    }
    
    label_position = 1.2*max(-log10(dataframe$pvalue))

    if (!is.na(subsample_point_cutoff)){
        above_cut = dataframe[-log10(pvalue) > subsample_point_cutoff]
        below_cut = dataframe[-log10(pvalue) <= subsample_point_cutoff] 
        set.seed(42)
        count = nrow(below_cut)*0.1
        rows = sample(nrow(below_cut), size = count)
        below_sample = below_cut[rows,]
        dataframe = rbind(above_cut, below_sample)
    }
    
    # scatter rows in dataframe to get scattered colors
    if (!is.na(gene_usage)){
        order = c('V-gene usage', 'J-gene usage', 'D-gene usage')
        dataframe$feature = factor(dataframe$feature, levels = order)
        dataframe = dataframe[order(dataframe$feature),]
    } else{
        set.seed(42)
        rows = sample(nrow(dataframe))
        dataframe = dataframe[rows,]
    }

    snps = ggplot(dataframe) +
        geom_vline(data = GENE_ANNOTATIONS, aes(xintercept = mid_pos), size = 1, lty = 'longdash', alpha = 0.4)+
        geom_text(data = GENE_ANNOTATIONS, size = 8, angle = 90, y = label_position, aes(x = mid_pos, label = gene_locus), vjust = 1.2, family = 'Arial') + 
        geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = feature), alpha = 0.75, size = point_size) +
        # geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.6) +
        facet_grid(productivity~chr, switch="x", space='free_x', scales = "free_x") +
        theme_classic() +
        theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.text.x = element_blank(), strip.text.x = element_text(size = 18)) +
        labs(y="-log10(p-value)", x="Chromosome Position") +
        geom_hline(data = sig_lines[1,], aes(yintercept=-log10(cutoff)), color = 'grey70', size=2.5) +
        ggtitle(plot_title) +
        scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8), expand = c(.4, .4)) +
        scale_color_manual(guide = guide_legend(reverse=TRUE), values = feature_colors) +
        ylim(c(min(-log10(dataframe$pvalue)), label_position*1.1)) + 
        scale_linetype_manual(name = 'significance', values = 2, guide = guide_legend(override.aes = list(color = 'grey70')))
    
    # ggsave(paste0(file_name, '.pdf'), plot = snps, width = 25, height = 18, units = 'in', dpi = 750, device = cairo_pdf)
    saveRDS(snps, file = paste0(file_name, '.rds'))
}
      
manhattan_plot_gene <- function(dataframe, phenotype_values = NA, plot_title, file_name, gene_subset, plot_zoom = NA, plotting_features = 'gene', unique_colors = FALSE, unique_line = FALSE){
    if (!is.na(phenotype_values)){
        dataframe$feature = mapvalues(dataframe$phenotype, from = unique(dataframe$phenotype), to = phenotype_values) 
    } else {
        dataframe$feature = dataframe$phenotype
    }

    # set phenotype colors
    dataframe = dataframe[order(-feature)]
    if (isTRUE(unique_colors)){
        feature_colors = brewer.pal(n = 5, name = "Set2")[(5-length(unique(dataframe$feature))+1):5]
    } else {
        feature_colors = brewer.pal(n = max(4, length(unique(dataframe$feature))), name = "Set2")
    }
    names(feature_colors) = unique(dataframe$feature)

    significance_cutoff = determine_significance_cutoff(0.05, type = gene_subset, dataframe)

    GENE_ANNOTATIONS$mid_pos = rowMeans(GENE_ANNOTATIONS[, c('pos1', 'pos2')])

    gene = GENE_ANNOTATIONS[gene_common_name == gene_subset]
    dataframe = dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]

    if (plotting_features == 'gene_features'){
        plot_features = GENE_FEATURES[gene_common_name == gene_subset][order(pos1)]
    } else {
        plot_features = gene
        plot_features$name = plot_features$gene_locus
    }

    if (!is.na(plot_zoom)){
        zoom_region = GENE_ANNOTATIONS[gene_common_name == plot_zoom]
        plot_features = plot_features[pos1 < zoom_region$pos1, pos1 := zoom_region$pos1]
        plot_features = plot_features[pos2 > zoom_region$pos2, pos2 := zoom_region$pos2]

        dataframe = dataframe[chr == zoom_region$chr & hg19_pos > zoom_region$pos1 & hg19_pos < zoom_region$pos2]
    }

    #Set plot parameters
    point_size = 5 
   
    if (!('productive' %in% colnames(dataframe))) {
        dataframe$productivity = 'all'
    } else {
        dataframe[productive == TRUE, productivity := 'productive']
        dataframe[productive == FALSE, productivity := 'non-productive']
    }

    feature_colors = feature_colors[names(feature_colors) %in% unique(dataframe$feature)]
 
    # scatter rows in dataframe to get scattered colors
    set.seed(42)
    rows = sample(nrow(dataframe))
    dataframe = dataframe[rows,]

    label_position = 1.2*max(-log10(dataframe$pvalue))

    snps = ggplot(dataframe) +
        geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = feature), alpha = 0.75, size = point_size) +
        geom_rect(data = plot_features, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = name), color = NA, alpha = 0.1) +
        geom_text(data = plot_features, size = 8, angle = 90, y = label_position, aes(x = mid_pos, label = gene_locus), family = 'Arial') + 
        facet_grid(productivity~chr, switch="x", space='free_x', scales = "free_x") +
        theme_classic() +
        theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.text.x = element_blank()) +
        labs(y="-log10(p-value)", x="Chromosome Position", fill = paste0(unique(plot_features$gene_locus), ' locus')) +
        ggtitle(plot_title) +
        scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8)) +
        scale_color_manual(guide = guide_legend(reverse=TRUE), values = feature_colors) +
        scale_fill_grey(guide = FALSE)
    
    if (isTRUE(unique_line)){
    snps = snps + geom_hline(yintercept=-log10(significance_cutoff), color = "black", size=2.5, linetype = 'dashed') 
    } else {
        snps = snps + geom_hline(yintercept=-log10(significance_cutoff), color = "grey70", size=2.5)
    }
    # ggsave(paste0(file_name, '.pdf'), plot = snps, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    saveRDS(snps, file = paste0(file_name, '.rds'))

}
 
manhattan_plot_loci <- function(dataframe, phenotype_values = NA, plot_title, file_name, gene_subset, significance_cutoff_type = 'genome-wide', plot_zoom = NA){
    gene_features = GENE_FEATURES[gene_common_name == gene_subset][order(pos1)]

    if (!is.na(phenotype_values)){
        dataframe$feature = mapvalues(dataframe$phenotype, from = unique(dataframe$phenotype), to = phenotype_values) 
    } else {
        dataframe$feature = dataframe$phenotype
    }

    significance_cutoff = determine_significance_cutoff(0.05, significance_cutoff_type, dataframe)

    gene = GENE_ANNOTATIONS[gene_common_name == gene_subset]
    if (significance_cutoff_type != 'genome-wide'){
        dataframe = dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
        plot_title = paste0(plot_title, '\nGene level Bonferroni significance cutoff')
    # } else if (significance_cutoff == gene_subset){
    #     dataframe = dataframe[hg19_pos < (gene$pos2) & hg19_pos > (gene$pos1) & chr == gene$chr]
    #     plot_title = paste0(plot_title, '\nGene level Bonferroni significance cutoff')
    } else {
        max_pos = max(as.numeric(gene_features$pos2))
        min_pos = min(as.numeric(gene_features$pos1))
        dataframe = dataframe[hg19_pos < max_pos & hg19_pos > min_pos & chr == gene$chr]
        plot_title = paste0(plot_title, '\nGenome-wide Bonferroni significance cutoff')
    }

    sig_snps = dataframe[pvalue < significance_cutoff]

    if (!('productive' %in% colnames(dataframe))) {
        dataframe$productive = 'all'
    }
    sig_snps$productivity = ifelse(sig_snps$productive == FALSE, 'non-productive', 'productive')
    
    sig_snps$type = 'snp'
    gene_features$type = 'feature'
    if (significance_cutoff_type == gene_subset){
        gene_features$feature = paste0('extended\n', unique(gene_features$gene_locus), '\nlocus') 
    } else {
        gene_features$feature = paste0(unique(gene_features$gene_locus), '\nlocus')
    }
    gene_features$mid_pos = rowMeans(gene_features[, c('pos1', 'pos2')])
    
    if (!is.na(plot_zoom)){
        zoom_region = GENE_ANNOTATIONS[gene_common_name == plot_zoom]
        gene_features = gene_features[pos1 < zoom_region$pos1, pos1 := zoom_region$pos1]
        gene_features = gene_features[pos2 > zoom_region$pos2, pos2 := zoom_region$pos2]

        sig_snps = sig_snps[chr == zoom_region$chr & hg19_pos > zoom_region$pos1 & hg19_pos < zoom_region$pos2]
    }
   
    together = rbind(sig_snps, gene_features, fill = TRUE)
    together$feature = factor(together$feature, levels = unique(together$feature))

    feature_count = length(unique(together$name)[!is.na(unique(together$name))])
    feature_colors = brewer.pal(n = 8, name = "Set1")

    feature_colors = feature_colors[-2]

    together = together[order(-log10(pvalue))]
    
    plot = ggplot(together) + 
        facet_grid(type~., scales='free', space = 'free_y') +
        theme_minimal() +
        theme(text = element_text(size = 40), axis.text.x = element_blank(), axis.title.y = element_blank()) +
        labs(x="Position") +
        ggtitle(plot_title) +
        geom_point(data = subset(together, type == 'snp'), aes(x = hg19_pos, y = productivity, color = -log10(pvalue)), size = 12, alpha = 0.8) + 
        facet_wrap(~feature, ncol = 1, strip.position = 'left', scales = 'free_y') +
        scale_color_gradient(low = "white", high = 'royalblue4', guide = guide_colorbar(barheight = 15), limits = c(significance_cutoff, max(-log10(together$pvalue)))) +
        new_scale_color() +
        geom_point(data = subset(together, type == 'feature'), aes(x = pos1, y = name, col = name), size = 4, alpha = 0.7, shape = '\u25BC', stroke = 5, position = position_nudge(y = 0.08)) +
        geom_point(data = subset(together, type == 'feature'), aes(x = pos2, y = name, col = name), size = 4, alpha = 0.7, shape = '\u25B2', stroke = 5, position = position_nudge(y = -0.08)) +
        geom_segment(data = subset(together, type == 'feature'), aes(x = pos1, xend = pos2, y = name, yend = name, color = name), alpha = 0.7, size = 3) +
        scale_color_manual(guide = 'none', values = feature_colors)
        
    final = plot + force_panelsizes(rows = c(1, 1, 1))

    # ggsave(paste0(file_name, '.pdf'), plot = final, width = 35, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    saveRDS(final, file = paste0(file_name, '.rds'))
}
    
set_plot_title_compare <- function(){ 
    if (PHENOTYPE == 'd0_trim'){
        title = paste0(PHENOTYPE_VALUE, ' trimming')
    }
    return(title)
}

get_file_name_compare <- function(gene){
    file = paste0(PROJECT_PATH, '/tcr-gwas/figures/', PHENOTYPE_TYPE, '_', gene, '_compare_with_without_allele_linkage_', CORRECTION_TYPE, 'correction')
    return(file)
}

manhattan_plot_gene_compare <- function(dataframe_correction, dataframe_no_correction, phenotype_values = NA, plot_title, file_name, gene_subset, plot_zoom = NA, plotting_features = 'gene'){
    dataframe_correction$linkage_correction = 'With correction'
    dataframe_no_correction$linkage_correction = 'No correction'
    dataframe = rbind(dataframe_correction, dataframe_no_correction)
    significance_cutoff = determine_significance_cutoff(0.05, type = gene_subset, dataframe_correction)

    if (!is.na(phenotype_values)){
        dataframe$feature = mapvalues(dataframe$phenotype, from = unique(dataframe$phenotype), to = phenotype_values) 
    } else {
        dataframe$feature = dataframe$phenotype
    }

    # set phenotype colors
    dataframe = dataframe[order(-feature)]
    feature_colors = brewer.pal(n = max(4, length(unique(dataframe$feature))), name = "Set2")
    names(feature_colors) = unique(dataframe$feature)


    gene = GENE_ANNOTATIONS[gene_common_name == gene_subset]
    dataframe = dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]

    if (plotting_features == 'gene_features'){
        plot_features = GENE_FEATURES[gene_common_name == gene_subset][order(pos1)]
    } else {
        plot_features = gene
        plot_features$name = plot_features$gene_locus
    }

    if (!is.na(plot_zoom)){
        zoom_region = GENE_ANNOTATIONS[gene_common_name == plot_zoom]
        plot_features = plot_features[pos1 < zoom_region$pos1, pos1 := zoom_region$pos1]
        plot_features = plot_features[pos2 > zoom_region$pos2, pos2 := zoom_region$pos2]

        dataframe = dataframe[chr == zoom_region$chr & hg19_pos > zoom_region$pos1 & hg19_pos < zoom_region$pos2]
    }

    #Set plot parameters
    point_size = 5 
   
    if (!('productive' %in% colnames(dataframe))) {
        dataframe$productivity = 'all'
    } else {
        dataframe[productive == TRUE, productivity := 'productive']
        dataframe[productive == FALSE, productivity := 'non-productive']
    }

    feature_colors = feature_colors[names(feature_colors) %in% unique(dataframe$feature)]
 
    # scatter rows in dataframe to get scattered colors
    set.seed(42)
    rows = sample(nrow(dataframe))
    dataframe = dataframe[rows,]
 
    snps = ggplot(dataframe) +
        geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = feature), alpha = 0.75, size = point_size) +
        geom_rect(data = plot_features, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = name), alpha = 0.1) +
        facet_grid(productivity ~ linkage_correction)+
        theme_classic() +
        theme(text = element_text(size = 40), axis.text.x = element_blank()) +
        labs(y="-log10(p-value)", x="Chromosome Position", fill = paste0(unique(plot_features$gene_locus), ' locus')) +
        geom_hline(yintercept=-log10(significance_cutoff), color = "grey70", size=2) +
        ggtitle(plot_title) +
        scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8)) +
        scale_color_manual(guide = guide_legend(reverse=TRUE), values = feature_colors) +
        scale_fill_grey()

    # ggsave(paste0(file_name, '.pdf'), plot = snps, width = 18, height = 20, units = 'in', dpi = 750, device = cairo_pdf)
    saveRDS(snps, file = paste0(file_name, '.rds'))
}
 
