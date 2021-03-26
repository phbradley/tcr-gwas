map_scanID_to_localID <- function(scanIDs_to_convert){
    ID_map_file = fread(ID_MAPPING_FILE)
    converted_IDs = plyr::mapvalues(scanIDs_to_convert, 
                                    ID_map_file$scanID, 
                                    ID_map_file$localID)
    return(converted_IDs)
}

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][8]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}



compile_all_genotypes_snp_list <- function(snp_list){
    require(SNPRelate)
    snp_gds_file = openfn.gds(SNP_GDS_FILE)
    genotypes = snpgdsGetGeno(snp_gds_file, snp.id = snp_list, with.id = TRUE)
    closefn.gds(snp_gds_file)
    
    genotype_matrix = genotypes$genotype
    rownames(genotype_matrix) = map_scanID_to_localID(genotypes$sample.id)
    colnames(genotype_matrix) = genotypes$snp.id 
    genotype_dt = data.table(localID = row.names(genotype_matrix), genotype_matrix)
    return(genotype_dt)
}

get_d_gene_usage <- function(){
    require(doParallel)
    require(foreach)
    files = list.files(TCR_REPERTOIRE_DATA_DIRECTORY, pattern = "*.tsv", full.names=TRUE)
    d_usage_data = data.table()        

    registerDoParallel(cores=NCPU)
    d_usage_data = foreach(file = files, .combine = 'rbind') %dopar% {
        file_data = fread(file)
        counts = file_data[,.N, by = d_gene]
        counts = counts[, proportion := N/sum(N)]
        counts$localID = extract_subject_ID(file)
        counts
    }
    stopImplicitCluster()
    return(d_usage_data) 
}

set_plot_title <- function(){
    return(paste0('SNP genotype versus TRBD2 alternate allele usage'))
}

set_file_name <- function(snp){
    path = paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/d_gene_allele_tcrb_analysis/')
    name = paste0('TRBD2_alt_allele_usage_by_SNP_genotype_', substring(snp, 4), '.pdf')
    return(paste0(path, name))
}

set_plot_subtitle <- function(snpID, sig_snps){
    snp_rows = sig_snps[snp == as.numeric(substring(snpID, 4))]
    snp_rows$productivity = ifelse(snp_rows$productive == 'TRUE', 'productive TCRs', 'non-productive TCRs')
    subtitle = paste0(snp_rows$phenotype, ' of ', snp_rows$productivity, ' pvalue: ', signif(snp_rows$pvalue, digits = 3))
    if (length(subtitle) > 1){
        paste(subtitle, collapse = '\n')
    }
    return(subtitle)
}

boxplot_by_d_allele_usage <- function(snp, dataframe, sig_snps){
    require(ggplot2)
    require(ggpubr)
    title = set_plot_title()
    subtitle = set_plot_subtitle(snp, sig_snps)
    filename = set_file_name(snp)

    compare = list(c('0', '1'), c('1', '2'), c('0', '2'))
    plot = ggboxplot(dataframe, x = snp, y = 'alt_allele_prop', fill = snp, size = 1.5, outlier.shape = NA) +
         stat_compare_means(comparisons = compare, size = 8) +
         # stat_compare_means(label.y = 0.85, size = 8) +
         geom_jitter(shape=16, position=position_jitter(0.05), size = 4, alpha = 0.75) +
         theme_classic() +
         theme(text = element_text(size = 30), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
         labs(title = title, subtitle = subtitle, y = 'TRBD2*02 allele usage')

     ggsave(filename, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = 'pdf')
}

######################################################
## functions to determine significance snps by gene ##
######################################################
get_gene_region_associations <- function(dataframe, gene, type = 'gene'){
    gene_dt = GENE_ANNOTATIONS[gene_common_name == gene]
    if (type == 'gene'){
        dataframe = dataframe[hg19_pos < (gene_dt$pos2) & hg19_pos > (gene_dt$pos1) & chr == gene_dt$chr]
    }
    else {
        dataframe = dataframe[hg19_pos < (gene_dt$pos2 + 200000) & hg19_pos > (gene_dt$pos1 - 200000) & chr == gene_dt$chr]
    }
    return(dataframe)
}

determine_stat_significance_cutoff <- function(alpha, significance_cutoff_type, genome_wide_dataframe, gene = NA, phenotype = unique(genome_wide_dataframe$phenotype)){
    stopifnot(significance_cutoff_type %in% c('genome-wide', 'gene', 'gene-surround'))
    features = length(phenotype)*length(unique(genome_wide_dataframe$productive))
    if (significance_cutoff_type == 'genome-wide'){
        snps_total = length(unique(genome_wide_dataframe$snp))
    } else {
        stopifnot(!is.na(gene))
        gene_dt = GENE_ANNOTATIONS[gene_common_name == gene]
        if (significance_cutoff_type == 'gene'){
            dataframe = genome_wide_dataframe[hg19_pos < (gene_dt$pos2) & hg19_pos > (gene_dt$pos1) & chr == gene_dt$chr]
        }
        else {
            dataframe = genome_wide_dataframe[hg19_pos < (gene_dt$pos2 + 200000) & hg19_pos > (gene_dt$pos1 - 200000) & chr == gene_dt$chr]
        }

        snps_total = length(unique(dataframe$snp))
    }
    sig_cutoff = (alpha/snps_total)/features
    return(sig_cutoff)
}


look_for_feature_overlap <- function(dataframe){
    for (feature in nrow(GENE_FEATURES)){
        gene_feature = GENE_FEATURES[feature,]
        dataframe[hg19_pos < gene_feature$pos2 & hg19_pos > gene_feature$pos1 & chr == gene_feature$chr, feature:= paste0(gene_feature$gene_locus, ': ', gene_feature$name)]
    }
    return(dataframe)
}

calculate_lambda_by_phenotype <- function(phenotype_dataframe){
    phenotype_dataframe$chisq = qchisq(1-phenotype_dataframe$pvalue, 1)
    lambda = median(phenotype_dataframe$chisq)/qchisq(0.5,1)
    return(lambda)
}

get_sig_snps_stats <- function(dataframe, name, gene, significance_cutoff_type_for_gene_locus){
    if (name == 'gene_usage'){
        significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = 'genome-wide', dataframe, phenotype = unique(dataframe$gene))
    } else {
        significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = 'genome-wide', dataframe)
    }
    print(paste0('The significance cutoff for ', name, ' is ', significance_cutoff))
    
    print(paste0('There are ', nrow(dataframe[pvalue < significance_cutoff]), ' significant associations for ', name, ' at a genome-wide significance threshold'))
    print(head(dataframe[pvalue < significance_cutoff][order(pvalue)]))

    gene_dataframe = get_gene_region_associations(dataframe, gene, type = significance_cutoff_type_for_gene_locus)
    if (name == 'gene_usage'){
        gene_significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = significance_cutoff_type_for_gene_locus, dataframe, gene = gene, phenotype = unique(dataframe$gene))
    } else {
        gene_significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = significance_cutoff_type_for_gene_locus, dataframe, gene = gene)
    }

    # gene_dataframe = look_for_feature_overlap(gene_dataframe)
    print(paste0('The significance cutoff for ', name, ' for ', gene, ' ', significance_cutoff_type_for_gene_locus, ' is ', gene_significance_cutoff))

    print(paste0('There are ', nrow(gene_dataframe[pvalue < gene_significance_cutoff]), ' significant associations for ', name, ' at a ', gene, ' ', significance_cutoff_type_for_gene_locus, ' significance threshold'))

    print(gene_dataframe[pvalue < gene_significance_cutoff][order(pvalue)][1:10])
    print(gene_dataframe[pvalue < gene_significance_cutoff][, .N, by = .(phenotype, productive)])
    return(gene_dataframe[pvalue < gene_significance_cutoff])
}

process_names <- function(lambda_dataframe){
    short_names = c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'v_pnucs_fraction_zero_trimming_subset', 'd0_pnucs_fraction_zero_trimming_subset', 'd1_pnucs_fraction_zero_trimming_subset', 'j_pnucs_fraction_zero_trimming_subset', 'v_trim_naive', 'd0_trim_naive', 'd1_trim_naive', 'j_trim_naive', 'v_pnucs_count', 'd0_pnucs_count', 'd1_pnucs_count', 'j_pnucs_count')
    long_names = c('V-gene trimming', '5\'-end D-gene trimming', '3\'-end D-gene trimming', 'J-gene trimming', 'V-D-gene insertions', 'D-J-gene insertions', 'Fraction of untrimmed TCRs containing V-gene P-nucleotides', 'Fraction of untrimmed TCRs containing 5\'-end D-gene P-nucleotides', 'Fraction of untrimmed TCRs containing 3\'-end D-gene P-nucleotides', 'Fraction of untrimmed TCRs containing J-gene P-nucleotides', 'V-gene trimming, no gene conditioning', '5\'-end D-gene trimming, no gene conditioning', '3\'-end D-gene trimming, no gene conditioning', 'J-gene trimming, no gene conditioning', 'V-gene P-nucleotide count', '5\'-end D-gene P-nucleotide count', '3\'-end D-gene P-nucleotide count', 'J-gene P-nucleotide count')
    names(long_names) = short_names

    long_names_present = long_names[names(long_names) %in% unique(lambda_dataframe$feature)]

    lambda_dataframe$feature = mapvalues(lambda_dataframe$feature, from = names(long_names_present), to = long_names_present) 

    return(lambda_dataframe)
}

get_lambdas_xtable <- function(dataframe, name){
    together = data.table()
    for (feature in unique(dataframe$phenotype)){
        for (prod in unique(dataframe$productive)){
            prod_long = ifelse(prod == TRUE, 'productive', 'non-productive')
            lambda = calculate_lambda_by_phenotype(dataframe[phenotype == feature & productive == prod])
            together = rbind(together, data.table(feature = feature, productivity = prod_long, lambda = lambda))
            print(paste0('The genome-wide lambda value for ', feature, ' and productivity=', prod, ' is ', lambda))
        }
    }
    together = process_names(together)
    latex_table = xtable(together)
    print(latex_table, file = paste0('analysis/', name, '_lambda_latex_table.txt'), include.rownames = FALSE)
    return(together)
}

get_lambdas_xtable_gene_usage <- function(dataframe){
    together = data.table()
    for (gene_name in unique(dataframe$gene)){
        for (prod in unique(dataframe$productive)){
            prod_long = ifelse(prod == TRUE, 'productive', 'non-productive')
            lambda = calculate_lambda_by_phenotype(dataframe[gene == gene_name & productive == prod])
            together = rbind(together, data.table(feature = gene_name, productivity = prod_long, lambda = lambda))
            print(paste0('The genome-wide lambda value for ', gene_name, ' and productivity=', prod, ' is ', lambda))
        }
    }
    latex_table = xtable(together)
    print(latex_table, file = paste0('analysis/gene_usage_lambda_latex_table.txt'))
    return(together)
}



######################################################
### D gene allele SNP linkage correction file and plotting ###
######################################################

get_allele_linkage_file_name <- function(method, gene){
    file_name = paste0(OUTPUT_PATH, '/results/d_allele_linkage/', PHENOTYPE, '_regressions_', CONDENSING_VARIABLE, '_d_infer-', INFER_MISSING_D_GENE, '_', PCA_COUNT, '_PCAir_PCs_', gene, '_d_gene_allele_test_', method, '.tsv')
    return(file_name)
} 

set_plot_title_compare <- function(){
    if (PHENOTYPE == 'd0_trim'){
        title = paste0(PHENOTYPE_VALUE, ' trimming')
    }
    return(title)
}

get_file_name_compare <- function(gene){
    file = paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/', PHENOTYPE, '_', gene, '_compare_with_without_allele_linkage_correction.pdf')
    return(file)
}

manhattan_plot_gene_compare <- function(dataframe_correction, dataframe_no_correction, phenotype_values = NA, plot_title, file_name, gene_subset, plot_zoom = NA, plotting_features = 'gene'){
    dataframe_correction$linkage_correction = 'With correction'
    dataframe_no_correction$linkage_correction = 'No correction'
    dataframe = rbind(dataframe_correction, dataframe_no_correction)

    if (!is.na(phenotype_values)){
        dataframe$feature = mapvalues(dataframe$phenotype, from = unique(dataframe$phenotype), to = phenotype_values) 
    } else {
        dataframe$feature = dataframe$phenotype
    }

    # set phenotype colors
    dataframe = dataframe[order(-feature)]
    feature_colors = brewer.pal(n = max(4, length(unique(dataframe$feature))), name = "Set2")
    names(feature_colors) = unique(dataframe$feature)

    significance_cutoff = determine_significance_cutoff(0.05, type = gene_subset, dataframe)

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
        geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = feature, shape = productivity), alpha = 0.5, size = point_size) +
        geom_rect(data = plot_features, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = name), alpha = 0.1) +
        facet_grid(.~linkage_correction)+
        theme_classic() +
        theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 40), axis.text.x = element_blank()) +
        labs(y="-log10(p-value)", x="Chromosome Position", fill = paste0(unique(plot_features$gene_locus), ' locus')) +
        geom_hline(yintercept=-log10(significance_cutoff), color = "grey70", size=2) +
        ggtitle(plot_title) +
        scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8)) +
        scale_color_manual(guide = guide_legend(reverse=TRUE), values = feature_colors) +
        scale_fill_grey()

    ggsave(file_name, plot = snps, width = 18, height = 10, units = 'in', dpi = 750, device = 'pdf')
}
 
