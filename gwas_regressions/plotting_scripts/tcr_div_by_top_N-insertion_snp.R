library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(SNPRelate)

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

tcr_div = fread(file = paste0(PROJECT_PATH, '/tcr-gwas/_ignore/emerson_tcrdiv_scores.tsv'))
setnames(tcr_div, 'hipid', 'localID')
setnames(tcr_div, 'tcrdiv', 'tcr_div')

top_insertion_snps = compile_manhattan_plot_data(c('vd_insert', 'dj_insert', 'total_insert'))[pvalue < 0.0005]

most_sig_snp = top_insertion_snps[order(pvalue)][1]

PHENOTYPE = most_sig_snp$phenotype
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/regression_functions.R'))

snp_gds_file = snpgdsOpen(SNP_GDS_FILE)
genotypes = as.data.table(snpgdsGetGeno(snp_gds_file, snp.id = most_sig_snp$snp, with.id = TRUE))
closefn.gds(snp_gds_file)
genotypes$localID = map_scanID_to_localID(genotypes$sample.id)
genotypes$genotype = as.character(genotypes$genotype)
together = merge(tcr_div, genotypes, by = 'localID')
together = together[!is.na(genotype)]
# average_all = data.table(tcr_div_mean = mean(together$tcr_div))
genotype_levels = list( c("0", "1"), c("1", "2"), c("0", "2") )
ggplot(together, aes(x=genotype, y=tcr_div, fill = genotype)) +
    geom_boxplot(lwd = 1.5) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 4, alpha = 0.5) +
    # stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", size = 10) +
    stat_compare_means(comparisons = genotype_levels, size = 10) +
    theme_classic() + 
    theme(text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
    ggtitle(paste0('TCR Repertoire Diversty by Allele Count')) +  
    # geom_hline(data = average_all, aes(yintercept = tcr_div_mean), size = 2.5, color = 'red', linetype = 2) +
    xlab(paste0('Allele Count for SNP ', most_sig_snp$snp)) +
    ylab('TCR Diversity') 

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/tcr_div_by_most_sig_N-insertion_allele.pdf'), plot = last_plot(), width = 15, height = 12, units = 'in', dpi = 750, device = 'pdf')



