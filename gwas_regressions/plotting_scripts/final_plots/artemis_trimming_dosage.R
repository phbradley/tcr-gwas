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

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))

# %PB not sure where to put this, but it might be nice to add a supplementary
# figure showing the trimming distributions for a common gene as a function of
# Artemis SNP genotype, something like the one I posted t    o slack. Where you
# can see the ``dosage'' effect of the minor allele. And also, truth in
# advertising, the relatively small magnitude of the overall change.
create_distribution_data <- function(data, feature_of_interest,gene_of_interest){
    gene_class = paste0(tolower(substring(gene_of_interest, 4, 4)), '_gene')
    subset = data[get(gene_class) == gene_of_interest, mean(get(feature_of_interest)), by = .(productive, localID)]
    setnames(subset, 'V1', 'mean_trim')
    return(subset)
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

boxplot_by_snp <- function(trimming_data, genotype_data, gene_of_interest, feature_of_interest){
    mean_trimming_by_gene = create_distribution_data(trimming_data, feature_of_interest, gene_of_interest)
    together = merge(mean_trimming_by_gene, genotype_data, by = 'localID')
    
    filtered = together[!is.na(snp)]
    filtered$snp = factor(filtered$snp, levels = c('0', '1', '2')) 
    filtered[productive == TRUE, productivity := 'productive']
    filtered[productive == FALSE, productivity := 'non-productive']

    print(cor.test(as.numeric(filtered[productive == TRUE]$snp), filtered[productive == TRUE]$mean_trim, method = "pearson")) 
    print(cor.test(as.numeric(filtered[productive == FALSE]$snp), filtered[productive == FALSE]$mean_trim, method = "pearson"))

    plot = ggboxplot(filtered, x = 'snp', y = 'mean_trim', fill = 'snp', lwd = 1.5) +
        facet_wrap(~productivity)+
        geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
        # stat_compare_means(comparisons = comparisons, aes(label = paste0('p = ', ..p.format..)), size = 10, family = 'Arial') +
        # stat_pvalue_manual(t_test, label = 'p', size = 10, family = 'Arial') +
        theme_classic(base_family = 'Arial') + 
        theme(text = element_text(size = 40, family = 'Arial'),legend.position = "none") +
        scale_fill_brewer(palette="Greys")

    final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 35), axis.text.x=element_text(size = 16), axis.text.y = element_text(size = 16), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'y') + panel_border(color = 'gray60', size = 1.5) 
    return(final_plot)
}


trimmings = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/trim_by_patient.tsv')[,-1]
setnames(trimmings, 'patient_id', 'localID')

trimming_associations = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))
gene = GENE_ANNOTATIONS[gene_common_name == 'artemis']

artemis_associations = trimming_associations[hg19_pos < (gene$pos2 + 200000) & hg19_pos > (gene$pos1 - 200000) & chr == gene$chr]
 
top_associations = artemis_associations[order(pvalue)][1:10]
top_associations[, min_p := min(pvalue), by = .(phenotype, productive)]

top_j_snp = top_associations[phenotype == 'j_trim'][1]
top_v_snp = top_associations[phenotype == 'v_trim'][1]

j_genotypes = compile_all_genotypes_snp_list(top_j_snp$snp)
v_genotypes = compile_all_genotypes_snp_list(top_v_snp$snp)
colnames(j_genotypes)[-1] = paste0('snp', colnames(j_genotypes)[-1])
colnames(v_genotypes)[-1] = paste0('snp', colnames(v_genotypes)[-1])

j_genotypes = association_genotype_assignment(top_j_snp$slope, j_genotypes)
v_genotypes = association_genotype_assignment(top_v_snp$slope, v_genotypes)

v_gene_usage = trimmings[, .N, by = .(v_gene, localID, productive)]
v_gene_usage[, max_N := max(N), by = .(localID, productive)]

top_v_gene = v_gene_usage[max_N == N][, .N, by = .(v_gene)][order(-N)][1]

j_gene_usage = trimmings[, .N, by = .(j_gene, localID, productive)]
j_gene_usage[, max_N := max(N), by = .(localID, productive)]

#use the top TRBJ1 gene (so that there is no potential for D-gene misidentification biases)
top_j_gene = j_gene_usage[max_N == N][, .N, by = .(j_gene)][order(-N)][2]


j_trim_boxplot = boxplot_by_snp(trimmings, j_genotypes, top_j_gene$j_gene, 'j_trim')
final_j_trim_boxplot = j_trim_boxplot + xlab('Genotype of the top J-gene trimming associated DCLRE1C SNP') + ylab('Number of J-gene nucleotides deleted') + geom_smooth(formula = y~x, method = lm, aes(x=snp, y = mean_trim,group = productivity), size = 3, fill = '#4292c6', color = '#08519c')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/j_trim_boxplot_', top_j_gene$j_gene, '_snp', top_j_snp$snp, '.pdf'), plot = final_j_trim_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

v_trim_boxplot = boxplot_by_snp(trimmings, v_genotypes, top_v_gene$v_gene, 'v_trim')
final_v_trim_boxplot = v_trim_boxplot + xlab('Genotype of the top V-gene trimming associated DCLRE1C SNP') + ylab('Number of V-gene nucleotides deleted')+ geom_smooth(formula = y~x, method = lm, aes(x=snp, y = mean_trim,group = productivity), size = 3, fill = '#4292c6', color = '#08519c')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/v_trim_boxplot_', top_v_gene$v_gene, '_snp', top_v_snp$snp, '.pdf'), plot = final_v_trim_boxplot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

