library(cowplot)
library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(Cairo)

args = commandArgs(trailingOnly=TRUE)

NCPU <<- as.numeric(args[1])

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/maf_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

dntt = find_snp_start_by_position(chromosome = 10, position1 = 98064085, position2 =98098321)

ethnicity = fread(file = ETHNICITY)

snp_start = dntt[1]
count = dntt[2]

genotypes = compile_all_genotypes(as.numeric(snp_start), as.numeric(count))

genotype_dt = merge(ethnicity[,c('localID', 'race.g')], genotypes, by = 'localID')

insertions = compile_mean_phenotype_data(c('v_gene', 'd_gene', 'd_gene', 'j_gene'), c('vd_insert', 'vd_insert', 'dj_insert', 'dj_insert'))[,c('localID', 'vj_insert', 'vd_insert', 'dj_insert', 'productive')]

insertions = insertions[productive == TRUE,lapply(.SD, mean), by = .(localID, productive)]
together = merge(insertions, genotype_dt, by = 'localID')
together$total_inserts = together$vj_insert + together$vd_insert + together$dj_insert

maf_dt = data.table()
for (snp in colnames(genotypes)[colnames(genotypes) != 'localID']){
    minor_allele = determine_true_minor_allele(snp, together)
    for (race in c(unique(together$race.g), 'all')){
        maf = calculate_maf(snp, minor_allele, race, together)
        temp = data.table(snp = snp, ancestry_group = race, maf = maf, minor_allele = minor_allele)
        maf_dt = rbind(maf_dt, temp)
    }
}

maf_dt$ancestry_group = str_replace(maf_dt$ancestry_group, ' ', '\n')
maf_dt$pca_cluster = paste0('\"', maf_dt$ancestry_group, '\"-\nassociated')

bonferroni = 7.74e-9
sig_snps = find_significant_snps(c('vd_insert_no_pca', 'dj_insert_no_pca'), bonferroni)

sig_snps_together_mafs = maf_dt[snp %in% sig_snps$snp]
sig_snps_together_mafs[pca_cluster == '\"all\"-\nassociated', pca_cluster := 'population']

t_test = sig_snps_together_mafs %>%
    t_test(maf ~ pca_cluster, ref.group = 'population')

average_all = mean(sig_snps_together_mafs[pca_cluster == 'population']$maf)

plot = ggboxplot(sig_snps_together_mafs[pca_cluster != 'population'], x = 'pca_cluster', y = 'maf', fill = 'pca_cluster', lwd = 1.5) +
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
    stat_pvalue_manual(t_test, size = 8, y.position = 0.75, remove.bracket = TRUE, family = 'Arial', label = "p = {p}") +
    theme_classic() + 
    theme(text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
    ggtitle('Minor Allele Frequency by Ancestry Group') +
    geom_hline(yintercept = average_all, size = 2.5, color = 'red', linetype = 2) +
    xlab('PCA cluster') +
    ylab('Minor allele frequency') + 
    scale_fill_brewer(palette = 'Set2')

final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 35), axis.text.x=element_text(size = 22), axis.text.y = element_text(size = 22), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'xy') + ylim(-0.01,0.85) + panel_border(color = 'gray60', size = 1.5)

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/figures/maf_by_race_boxplot.pdf'), plot = final_plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

