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

get_actual_sample_size_by_genotypes_snp_list <- function(snp_list){
    genotypes = compile_all_genotypes_snp_list(snp_list)
    sample_size = as.data.frame(colSums(!is.na(genotypes)))
    sample_size$snp = rownames(sample_size)
    rownames(sample_size) = NULL
    colnames(sample_size) = c('sample_size', 'snp')
    return(sample_size)
}

snp_file_by_snp_list <- function(snp_list){
    snp_data = fread(SNP_META_DATA_FILE)[, -c('V1')]
    snp_data_subset = snp_data[snpid %in% snp_list]
    colnames(snp_data_subset) = c('snp', 'hg19_pos', 'chr', 'snpallele')
    return(snp_data_subset[, -c('snpallele')])
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
    return(paste0('SNP genotype versus TRBD2*02 usage'))
}

set_file_name <- function(snp, type = 'allele_usage', final_figure = FALSE){
    if (final_figure == TRUE){
        path = paste0(PROJECT_PATH, '/tcr-gwas/figures/')
    } else {
        path = paste0(PROJECT_PATH, '/tcr-gwas/figures/d_gene_allele_tcrb_analysis/')
    }
    if (type == 'allele_usage'){
        name = paste0('TRBD2_alt_allele_usage_by_SNP_genotype_', substring(snp, 4), '.pdf')
    } else if (type == 'allele_genotype'){
        name = paste0('TRBD2_allele_genotype_by_SNP_genotype_', substring(snp, 4), '.pdf')
    }
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

boxplot_by_d_allele_usage <- function(snp, dataframe, sig_snps, final_figure = FALSE){
    require(ggplot2)
    require(ggpubr)
    title = set_plot_title()
    subtitle = set_plot_subtitle(snp, sig_snps)
    filename = set_file_name(snp, type = 'allele_usage', final_figure)
    dataframe = dataframe[!is.na(get(snp))]
    compare = list(c('0', '1'), c('1', '2'), c('0', '2'))
    plot = ggboxplot(dataframe, x = snp, y = 'alt_allele_prop', fill = snp, size = 2, outlier.shape = NA) +
         # stat_compare_means(label.y = 0.85, size = 8) +
         geom_jitter(shape=16, position=position_jitter(0.05), size = 6, alpha = 0.75) +
         theme_classic() +
         theme(text = element_text(size = 40), legend.position = "none") +
         scale_fill_brewer(palette = 'Set2')
    
    if (final_figure == TRUE){
        final_plot = plot + 
            # stat_compare_means(label='p.signif', comparisons = compare, size = 10, method = 't.test', p.adjust.method = "bonferroni", family = 'Arial') +
            labs(title = title, x = 'Top association TCRB SNP', y = 'TRBD2*02 allele usage')
    } else {
        final_plot = plot + 
         stat_compare_means(comparisons = compare, size = 10, method = 't.test', p.adjust.method = "bonferroni") +
         labs(title = title, subtitle = subtitle, y = 'TRBD2*02 allele usage')
    }
    return(final_plot)
     ggsave(filename, plot = final_plot, width = 14, height = 10, units = 'in', dpi = 750, device = 'pdf')
}

boxplot_by_d_allele_genotype <- function(snp, dataframe, sig_snps, final_figure = FALSE){
    require(ggplot2)
    require(ggpubr)
    title = set_plot_title()
    subtitle = set_plot_subtitle(snp, sig_snps)
    filename = set_file_name(snp, type = 'allele_genotype', final_figure)
    dataframe = dataframe[!is.na(get(snp))]
    compare = list(c('0', '1'), c('1', '2'), c('0', '2'))
    plot = ggboxplot(dataframe, x = snp, y = 'TRBD2_alt_allele_genotype', fill = snp, size = 2, outlier.shape = NA) +
         # stat_compare_means(label.y = 0.85, size = 8) +
         geom_point(shape=16, size = 6, alpha = 0.75, position = position_jitter(w = 0.1, h = 0)) +
         theme_classic() +
         theme(text = element_text(size = 40), legend.position = "none") +
         scale_fill_brewer(palette = 'Set2') +
         scale_y_continuous(breaks = c(0, 1, 2))
    
    if (final_figure == TRUE){
        final_plot = plot + 
            # stat_compare_means(label='p.signif', comparisons = compare, size = 10, method = 't.test', p.adjust.method = "bonferroni") +
            labs(title = title, x = 'Genotype of the top association TCRB SNP', y = 'TRBD2*02 allele genotype')
    } else {
        final_plot = plot + 
         # stat_compare_means(comparisons = compare, size = 10, method = 't.test', p.adjust.method = "bonferroni") +
         labs(title = title, subtitle = subtitle, y = 'TRBD2*02 allele genotype')
    }
    return(final_plot)
     ggsave(filename, plot = final_plot, width = 14, height = 10, units = 'in', dpi = 750, device = 'pdf')
}

heatmap_by_d_allele_genotype <- function(snp, dataframe, sig_snps, final_figure = FALSE){
    require(ggplot2)
    require(ggpubr)
    title = set_plot_title()
    subtitle = set_plot_subtitle(snp, sig_snps)
    filename = set_file_name(snp, type = 'allele_genotype', final_figure)
    dataframe = dataframe[!is.na(get(snp))]
    compare = list(c('0', '1'), c('1', '2'), c('0', '2'))
    plot = ggboxplot(dataframe, x = snp, y = 'TRBD2_alt_allele_genotype', fill = snp, size = 2, outlier.shape = NA) +
         # stat_compare_means(label.y = 0.85, size = 8) +
         geom_point(shape=16, size = 6, alpha = 0.75, position = position_jitter(w = 0.1, h = 0)) +
         theme_classic() +
         theme(text = element_text(size = 40), legend.position = "none") +
         scale_fill_brewer(palette = 'Set2') +
         scale_y_continuous(breaks = c(0, 1, 2))
    
    if (final_figure == TRUE){
        final_plot = plot + 
            # stat_compare_means(label='p.signif', comparisons = compare, size = 10, method = 't.test', p.adjust.method = "bonferroni") +
            labs(title = title, x = 'Genotype of the top association TCRB SNP', y = 'TRBD2*02 allele genotype')
    } else {
        final_plot = plot + 
         # stat_compare_means(comparisons = compare, size = 10, method = 't.test', p.adjust.method = "bonferroni") +
         labs(title = title, subtitle = subtitle, y = 'TRBD2*02 allele genotype')
    }
    return(final_plot)
     ggsave(filename, plot = final_plot, width = 14, height = 10, units = 'in', dpi = 750, device = 'pdf')
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

        gene_chr = gene_dt$chr
        gene_pos1 = gene_dt$pos1
        gene_pos2 = gene_dt$pos2

        if (significance_cutoff_type == 'gene'){
            dataframe = genome_wide_dataframe[hg19_pos < gene_pos2 & hg19_pos > gene_pos1 & chr == gene_chr]
        }
        else {
            dataframe = genome_wide_dataframe[hg19_pos < (gene_pos2 + 200000) & hg19_pos > (gene_pos1 - 200000) & chr == gene_chr]
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

get_scale_factor <- function(random_bootstrap_data, whole_genome_data, productivity){
    together = merge(whole_genome_data, random_bootstrap_data, by = c('snp', 'productive'))

    together = together[productive == productivity]
    together$z2.x = (together$slope.x/together$standard_error.x)^2
    together$z2.y = (together$slope.y/together$standard_error.y)^2

    scale = lm(together$z2.y ~ together$z2.x)
    return(scale)
}

temp_get_merged_data <- function(random_bootstrap_data, whole_genome_data, productivity){
    together = merge(whole_genome_data, random_bootstrap_data, by = c('snp', 'productive'))

    together = together[productive == productivity]
    together$z2.x = (together$slope.x/together$standard_error.x)^2
    together$z2.y = (together$slope.y/together$standard_error.y)^2

    return(together)
}


transform_z2_using_scale <- function(scale, z2){
    z2_transform = z2*(scale$coefficients[2]) + (scale$coefficients[1])
    return(z2_transform)
}


calculate_lambda_by_phenotype <- function(phenotype_dataframe){
    phenotype_dataframe$chisq = qchisq(1-phenotype_dataframe$pvalue, 1)
    lambda = median(phenotype_dataframe$chisq)/qchisq(0.5,1)
    return(lambda)
}

combine_rsids <- function(dataframe){
    rsids = fread(RSIDS)
    together = merge(dataframe, rsids)
    return(together)
}

get_snp_genotype <- function(rsid_name, phenotype){
    source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

    rsids = fread(RSIDS)
    rsid_info = rsids[substring(rsid, 1, nchar(rsid_name)) == rsid_name]
    phenotype_results = compile_manhattan_plot_data(phenotype)
    snp_id = rsid_info$snp
    phenotype_results_subset = phenotype_results[snp == snp_id]
    print(phenotype_results_subset)

    return(rsid_info)
}

get_association_count <- function(dataframe, name){
    if (name == 'gene_usage'){
        significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = 'genome-wide', dataframe, phenotype = unique(dataframe$gene))
    } else {
        significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = 'genome-wide', dataframe)
    }
    print(paste0('The significance cutoff for ', name, ' is ', significance_cutoff))
    
    print(paste0('There are ', nrow(dataframe[pvalue < significance_cutoff]), ' significant associations for ', name, ' at a genome-wide significance threshold'))

    sigs = dataframe[pvalue < significance_cutoff]
    sigs = combine_rsids(sigs)

    print(sigs[, .N, by = .(phenotype, productive)])

    return(sigs)
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

    gene_dataframe = combine_rsids(gene_dataframe)
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

get_lambdas_standard <- function(dataframe){
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
    return(together)
}

get_lambdas_gene_usage <- function(dataframe){
    together = data.table()
    for (gene_name in unique(dataframe$gene)){
        for (prod in unique(dataframe$productive)){
            prod_long = ifelse(prod == TRUE, 'productive', 'non-productive')
            lambda = calculate_lambda_by_phenotype(dataframe[gene == gene_name & productive == prod])
            together = rbind(together, data.table(feature = gene_name, productivity = prod_long, lambda = lambda))
            print(paste0('The genome-wide lambda value for ', gene_name, ' and productivity=', prod, ' is ', lambda))
        }
    }
    return(together)
}

get_lambdas <- function(dataframe, name){
    if (name == 'gene_usage'){
        lambdas = get_lambdas_gene_usage(dataframe)
    } else {
        lambdas = get_lambdas_standard(dataframe)
    }

    fwrite(lambdas, file = paste0(OUTPUT_PATH, '/results/lambdas/', name, '_lambda_table.txt'), sep = ',')
    return(lambdas)
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
    file = paste0(PROJECT_PATH, '/tcr-gwas/figures/', PHENOTYPE, '_', gene, '_compare_with_without_allele_linkage_', CORRECTION_TYPE, '_correction.pdf')
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

find_feature_overlaps <- function(dataframe, gene_subset, significance_cutoff_type = 'genome-wide'){
    gene_features = GENE_FEATURES[gene_common_name == gene_subset][order(pos1)]

    significance_cutoff = determine_significance_cutoff(0.05, significance_cutoff_type, dataframe)

    gene = GENE_ANNOTATIONS[gene_common_name == gene_subset]
    if (significance_cutoff_type != 'genome-wide'){
        dataframe = dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
    } else {
        max_pos = max(as.numeric(gene_features$pos2))
        min_pos = min(as.numeric(gene_features$pos1))
        dataframe = dataframe[hg19_pos < max_pos & hg19_pos > min_pos & chr == gene$chr]
    }

    sig_snps = dataframe[pvalue < significance_cutoff]
    print(paste0('total associations: ', nrow(sig_snps)))

    for (gene_feature in seq(1, nrow(gene_features))){
        feature = gene_features[gene_feature]
        chrom = as.numeric(feature$chr)
        pos1 = feature$pos1
        pos2 = feature$pos2
        if (gene_subset == 'tcrb'){
            name = feature$gene
        } else {
            name = feature$name
        }
        sig_snps[chr ==  chrom & hg19_pos > pos1 & hg19_pos < pos2, feature := name]
    }

    # sig_snps[chr == gene$chr | hg19_pos < pos1 | hg19_pos > pos2, feature := 'not in the locus']
    print(sig_snps[, .N, by = feature])
}

clean_supp_data <- function(dataframe){
    dataframe = combine_rsids(dataframe)
    sample_sizes = get_actual_sample_size_by_genotypes_snp_list(dataframe$snp)
    together = merge(dataframe, sample_sizes)
    final = together[, c('rsid', 'chr', 'hg19_pos', 'phenotype', 'productive', 'slope', 'standard_error', 'pvalue', 'sample_size')]
    colnames(final) = c('rsid', 'chr', 'hg19_pos', 'phenotype', 'productive', 'beta', 'standard_error', 'pvalue', 'sample_size')
    return(final)
}

compile_random_bootstrapped_results <- function(){
    path = paste0(OUTPUT_PATH, '/results/bootstrap_lambda_analysis/', PHENOTYPE, '_', REPETITIONS)
    files = fs::dir_ls(path)
    if (length(files) <= 1){
        stop('Need to run `submit_all_random_bootstraps_analyses.sh` for specified phenotype before calculating lambda')
    }
    compiled_filename = paste0(OUTPUT_PATH, '/results/bootstrap_lambda_analysis/', PHENOTYPE, '_', REPETITIONS, '_random_bootstraps.tsv')
    vroom::vroom_write(vroom::vroom(files, num_threads = NCPU), compiled_filename)
}

get_random_boostrap_results <- function(){
    filename = paste0(OUTPUT_PATH, '/results/bootstrap_lambda_analysis/', PHENOTYPE, '_', REPETITIONS, '_random_bootstraps.tsv')
    if (!file.exists(filename)){
        compile_random_bootstrapped_results()
    }
    file = fread(filename)
    return(file)
}
