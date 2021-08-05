compile_data <- function(phenotype){
    PHENOTYPE <<- phenotype
    source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/validation_cohort_analysis/regression_functions_validation_cohort.R'))
    phenotypes = compile_phenotype_data()
    genotypes = compile_all_genotypes(VALIDATION_SNP_PATH)
    colnames(genotypes) = c('localID', 'rs12768894', 'rs3762093')
    together = merge(phenotypes, genotypes, by = 'localID', all.x = TRUE)
    return(together)
}

order_allele_genotypes <- function(genotype_data){
    ordered_genotypes = genotype_data[, .N, by = .(snp, allele_genotype)][order(snp)]$allele_genotype
    genotype_data$allele_genotype = factor(genotype_data$allele_genotype, levels = ordered_genotypes)
    return(genotype_data)
}

association_genotype_assignment <- function(genotype_data, rs_id){
    alleles = fread(VALIDATION_SNP_ALLELES)
    colnames(alleles) = c('localID', 'rs12768894_alleles', 'rs3762093_alleles')
    genotype_data = merge(genotype_data, alleles)
    genotype_data$switched_genotype = FALSE
    genotype_data$allele_genotype = genotype_data[[paste0(rs_id, '_alleles')]] 
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

get_regression_data <- function(snp_name, phenotype){
    PHENOTYPE <<- phenotype
    source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/validation_cohort_analysis/regression_functions_validation_cohort.R'))
    file = make_regression_file_name()
    data = fread(file)
    data = data[snp == snp_name]
    return(data)
}

make_title_from_regression_data <- function(data){
    stopifnot(length(unique(data$snp)) == 1)
    title = paste0(unique(data$snp), ' and ', unique(data$phenotype), '\nnon-productive p = ', signif(data[productive == FALSE]$pvalue, 3), ', productive p = ', signif(data[productive == TRUE]$pvalue, 3))
    return(title)
}

boxplot_by_snp <- function(data, snp, gene_of_interest = NA, feature_of_interest, final_plot = FALSE){
    #TODO incorporate gene of interest thing
    if (final_plot == FALSE){
        regression_data = get_regression_data(snp_name = snp, phenotype = feature_of_interest)
        title = make_title_from_regression_data(regression_data)
    }else{
        title = ''
    }

    filtered = data[!is.na(snp)]
    filtered$snp = factor(filtered[[snp]], levels = c('0', '1', '2')) 
    filtered[productive == TRUE, productivity := 'productive']
    filtered[productive == FALSE, productivity := 'non-productive']
    filtered = association_genotype_assignment(filtered, snp)

    stats = get_significance(filtered, feature_of_interest) 

    plot = ggboxplot(filtered, x = 'allele_genotype', y = feature_of_interest, fill = 'snp', lwd = 1.5) +
        facet_wrap(~productivity)+
        geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
        # stat_compare_means(comparisons = comparisons, aes(label = paste0('p = ', ..p.format..)), size = 10, family = 'Arial') +
        # stat_pvalue_manual(t_test, label = 'p', size = 10, family = 'Arial') +
        # stat_pvalue_manual(stats, label = 'p', tip.length = 0.01, family = 'Arial', size = 8) +
        theme_classic(base_family = 'Arial') + 
        theme(text = element_text(size = 40, family = 'Arial'),legend.position = "none") +
        scale_fill_brewer(palette="Greys") +
        ggtitle(title) 

    final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 25), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + background_grid(major = 'y') + panel_border(color = 'gray60', size = 1.5) 
    return(final_plot)
}


