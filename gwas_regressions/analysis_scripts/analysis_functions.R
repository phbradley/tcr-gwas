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
    gene = GENE_ANNOTATIONS[gene_common_name == gene]
    if (type == 'gene'){
        dataframe = dataframe[hg19_pos < (gene$pos2) & hg19_pos > (gene$pos1) & chr == gene$chr]
    }
    else {
        dataframe = dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
    }
    return(dataframe)
}

determine_stat_significance_cutoff <- function(alpha, significance_cutoff_type, genome_wide_dataframe, gene = NA){
    stopifnot(significance_cutoff_type %in% c('genome-wide', 'gene', 'gene-surround'))
    features = length(unique(genome_wide_dataframe$phenotype))*length(unique(genome_wide_dataframe$productive))
    if (significance_cutoff_type == 'genome-wide'){
        snps_total = length(unique(genome_wide_dataframe$snp))
    } else {
        stopifnot(!is.na(gene))
        gene = GENE_ANNOTATIONS[gene_common_name == gene]
        if (significance_cutoff_type == 'gene'){
            dataframe = genome_wide_dataframe[hg19_pos < (gene$pos2) & hg19_pos > (gene$pos1) & chr == gene$chr]
        }
        else {
            dataframe = genome_wide_dataframe[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
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

get_sig_snps_stats <- function(dataframe, name, gene, significance_cutoff_type_for_gene_locus){
    significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = 'genome-wide', dataframe)
    print(paste0('The significance cutoff for ', name, ' is ', significance_cutoff))
    
    print(paste0('There are ', nrow(dataframe[pvalue < significance_cutoff]), ' significant associations for ', name, ' at a genome-wide significance threshold'))
    print(head(dataframe[pvalue < significance_cutoff][order(pvalue)]))

    gene_dataframe = get_gene_region_associations(dataframe, gene, type = significance_cutoff_type_for_gene_locus)
    gene_significance_cutoff = determine_stat_significance_cutoff(0.05, significance_cutoff_type = significance_cutoff_type_for_gene_locus, dataframe, gene = gene)
    
    # gene_dataframe = look_for_feature_overlap(gene_dataframe)
    print(paste0('The significance cutoff for ', name, ' for ', gene, ' ', significance_cutoff_type_for_gene_locus, ' is ', gene_significance_cutoff))

    print(paste0('There are ', nrow(gene_dataframe[pvalue < gene_significance_cutoff]), ' significant associations for ', name, ' at a ', gene, ' ', significance_cutoff_type_for_gene_locus, ' significance threshold'))

    print(gene_dataframe[pvalue < gene_significance_cutoff][order(pvalue)][1:10])
    print(gene_dataframe[pvalue < gene_significance_cutoff][, .N, by = .(phenotype, productive)])
    return(gene_dataframe[pvalue < gene_significance_cutoff])
}


