library(cowplot)
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(Cairo)

CHAIN <<- 'beta'
source('config/config.R')
source(paste0('config/validation_file_paths_', CHAIN, '.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))

compile_data <- function(phenotype){
    PHENOTYPE <<- phenotype
    source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions_validation_cohort.R'))
    phenotypes = compile_phenotype_data()
    genotypes = compile_all_genotypes(VALIDATION_SNP_PATH)
    colnames(genotypes) = c('localID', 'rs12768894', 'rs3762093')
    together = merge(phenotypes, genotypes, by = 'localID', all.x = TRUE)
    return(together)
}

association_genotype_assignment <- function(association_slope, genotype_data){
    snp_num = colnames(genotype_data)[colnames(genotype_data) != 'localID']
    if (association_slope < 0){
        genotype_data[get(snp_num)==0, snp := 2]
        genotype_data[get(snp_num)==2, snp := 0]
        genotype_data[get(snp_num)==1, snp := 1]
    }
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
    source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions_validation_cohort.R'))
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

    stats = get_significance(filtered, feature_of_interest) 

    plot = ggboxplot(filtered, x = 'snp', y = feature_of_interest, fill = 'snp', lwd = 1.5) +
        facet_wrap(~productivity)+
        geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
        # stat_compare_means(comparisons = comparisons, aes(label = paste0('p = ', ..p.format..)), size = 10, family = 'Arial') +
        # stat_pvalue_manual(t_test, label = 'p', size = 10, family = 'Arial') +
        stat_pvalue_manual(stats, label = 'p', tip.length = 0.01, family = 'Arial', size = 8) +
        theme_classic(base_family = 'Arial') + 
        theme(text = element_text(size = 40, family = 'Arial'),legend.position = "none") +
        scale_fill_brewer(palette="Greys") +
        ggtitle(title)

    final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + background_grid(major = 'y') + panel_border(color = 'gray60', size = 1.5) 
    return(final_plot)
}

vd_insertion_data = compile_data(phenotype = 'vd_insert')
dj_insertion_data = compile_data(phenotype = 'dj_insert')


# Make figures for insertions: 
vd_insert_boxplot = boxplot_by_snp(data = vd_insertion_data, snp = 'rs3762093', feature_of_interest = 'vd_insert')
final_vd_insert_boxplot = vd_insert_boxplot + xlab('rs3762093 SNP genotype') + ylab('Number of V-D-gene junction N-inserts')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/vd_insert_boxplot_rs3762093.pdf'), plot = final_vd_insert_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

dj_insert_boxplot = boxplot_by_snp(data = dj_insertion_data, snp = 'rs3762093', feature_of_interest = 'dj_insert')
final_dj_insert_boxplot = dj_insert_boxplot + xlab('rs3762093 SNP genotype') + ylab('Number of D-J-gene junction N-inserts')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/dj_insert_boxplot_rs3762093.pdf'), plot = final_dj_insert_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)


v_trim_data = compile_data(phenotype = 'v_trim')
j_trim_data = compile_data(phenotype = 'j_trim')
d0_trim_data = compile_data(phenotype = 'd0_trim')
d1_trim_data = compile_data(phenotype = 'd1_trim')

cdr3 = combine_genes_by_common_cdr3()
j_gene = 'TRBJ2-1*01'
j_gene_cdr3_id = cdr3[id == j_gene]$cdr3_gene_group
j_trim_data_subset = j_trim_data[cdr3_gene_group == j_gene_cdr3_id]

v_gene = 'TRBV12-3*01'
v_gene_cdr3_id = cdr3[id == v_gene]$cdr3_gene_group
v_trim_data_subset = v_trim_data[cdr3_gene_group == v_gene_cdr3_id]

d_gene = 'TRBD1*01'
d_gene_cdr3_id = cdr3[id == d_gene]$cdr3_gene_group
d0_trim_data_subset = d0_trim_data[cdr3_gene_group == d_gene_cdr3_id]
d1_trim_data_subset = d1_trim_data[cdr3_gene_group == d_gene_cdr3_id]

# Make figures for insertions: 
j_trim_boxplot = boxplot_by_snp(data = j_trim_data_subset, snp = 'rs12768894', feature_of_interest = 'j_trim')
final_j_trim_boxplot = j_trim_boxplot + xlab('rs12768894 SNP genotype') + ylab('Number of J-gene deletions')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/j_trim_boxplot_rs12768894.pdf'), plot = final_j_trim_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

v_trim_boxplot = boxplot_by_snp(data = v_trim_data_subset, snp = 'rs12768894', feature_of_interest = 'v_trim')
final_v_trim_boxplot = v_trim_boxplot + xlab('rs12768894 SNP genotype') + ylab('Number of J-gene deletions')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/v_trim_boxplot_rs12768894.pdf'), plot = final_v_trim_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

d0_trim_boxplot = boxplot_by_snp(data = d0_trim_data_subset, snp = 'rs12768894', feature_of_interest = 'd0_trim')
final_d0_trim_boxplot = d0_trim_boxplot + xlab('rs12768894 SNP genotype') + ylab('Number of J-gene deletions')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/d0_trim_boxplot_rs12768894.pdf'), plot = final_d0_trim_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

d1_trim_boxplot = boxplot_by_snp(data = d1_trim_data_subset, snp = 'rs12768894', feature_of_interest = 'd1_trim')
final_d1_trim_boxplot = d1_trim_boxplot + xlab('rs12768894 SNP genotype') + ylab('Number of J-gene deletions')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/d1_trim_boxplot_rs12768894.pdf'), plot = final_d1_trim_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)




