create_distribution_data <- function(data, feature_of_interest,gene_of_interest = NA){
    gene_class = paste0(tolower(substring(gene_of_interest, 4, 4)), '_gene')
    if (!is.na(gene_of_interest)){
        subset = data[get(gene_class) == gene_of_interest, mean(get(feature_of_interest)), by = .(productive, localID)]
    } else {
        subset = data[, mean(get(feature_of_interest)), by = .(productive, localID)]
    }
    setnames(subset, 'V1', 'mean')
    return(subset)
}

get_genotype_alleles <- function(snp_meta_data, switched_genotypes, genotypes){
    alleles = str_split(snp_meta_data$allele, '/')
    A_allele = alleles[[1]][1]
    B_allele = alleles[[1]][2]
    allele_genotypes = c()
    for (row in 1:length(genotypes)){
        if (is.na(genotypes[row])){
            allele_genotype = NA
        } else if (genotypes[row] == 0){
            allele_genotype = ifelse(switched_genotypes[row] == FALSE, paste0(B_allele, B_allele), paste0(A_allele, A_allele))
        } else if (genotypes[row] == 1){
            allele_genotype = ifelse(switched_genotypes[row] == FALSE, paste0(A_allele, B_allele), paste0(B_allele, A_allele))
        } else if (genotypes[row] == 2){
            allele_genotype = ifelse(switched_genotypes[row] == FALSE, paste0(A_allele, A_allele), paste0(B_allele, B_allele))
        }
        allele_genotypes = c(allele_genotypes, allele_genotype)
    }
    return(allele_genotypes)
}

order_allele_genotypes <- function(genotype_data){
    ordered_genotypes = genotype_data[, .N, by = .(snp, allele_genotype)][order(snp)]$allele_genotype
    genotype_data$allele_genotype = factor(genotype_data$allele_genotype, levels = ordered_genotypes)
    return(genotype_data)
}

association_genotype_assignment <- function(association_slope, genotype_data){
    snp_num = colnames(genotype_data)[colnames(genotype_data) != 'localID']
    snp_data = fread(SNP_META_DATA_FILE)
    snp_meta = snp_data[snpid == snp_num]

    if (association_slope < 0){
        genotype_data[get(snp_num)==0, snp := 2]
        genotype_data[get(snp_num)==2, snp := 0]
        genotype_data[get(snp_num)==1, snp := 1]
        genotype_data$switched_genotype = TRUE
    } else {
        genotype_data$snp = genotype_data[[snp_num]]
        genotype_data$switched_genotype = FALSE
    }

    genotype_data[, allele_genotype := get_genotype_alleles(snp_meta, switched_genotype, snp)]
    genotype_data = order_allele_genotypes(genotype_data)
    return(genotype_data)
}



get_significance <- function(data, feature_of_interest){
    stat.test = data %>% 
        group_by(productivity) %>%
        t_test(formula(paste0(feature_of_interest, '~ snp'))) %>% 
        add_significance() %>%
        add_xy_position(x = 'snp') %>%
        arrange(group2, n1, group1)
    return(stat.test)
}


boxplot_by_snp <- function(trimming_data, genotype_data, feature_of_interest, gene_of_interest = NA){
    mean_trimming_by_gene = create_distribution_data(trimming_data, feature_of_interest, gene_of_interest)
    together = merge(mean_trimming_by_gene, genotype_data, by = 'localID')
    
    filtered = together[!is.na(snp)]
    filtered$snp = factor(filtered$snp, levels = c('0', '1', '2')) 
    filtered[productive == TRUE, productivity := 'productive']
    filtered[productive == FALSE, productivity := 'non-productive']

    stats = get_significance(filtered, 'mean') 

    plot = ggboxplot(filtered, x = 'allele_genotype', y = 'mean', fill = 'snp', lwd = 1.5) +
        facet_wrap(~productivity)+
        geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
        # stat_pvalue_manual(stats, label = 'p', tip.length = 0.01, family = 'Arial', size = 8) +
        theme_classic(base_family = 'Arial') + 
        theme(text = element_text(size = 40, family = 'Arial'),legend.position = "none") +
        scale_fill_brewer(palette="Greys")

    final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 18), axis.text.y = element_text(size = 18), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'y') + panel_border(color = 'gray60', size = 1.5) 
    return(final_plot)
}


