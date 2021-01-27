library(rstatix)
library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(ggpubr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/maf_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

dntt = find_snp_start_by_position(chromosome = 10, position1 = 98064085, position2 =98098321)

ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))

snp_start = dntt[1]
count = dntt[2]

genotypes = compile_all_genotypes(as.numeric(snp_start), as.numeric(count))

genotype_dt = merge(ethnicity[,c('localID', 'race.g')], genotypes, by = 'localID')

insertions = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/insertions_by_patient.tsv')[,c('patient_id', 'vj_insert', 'vd_insert', 'dj_insert', 'productive')]

setnames(insertions, 'patient_id', 'localID')
insertions = insertions[productive == TRUE,lapply(.SD, mean), by = .(localID, productive)]
together = merge(insertions, genotype_dt, by = 'localID')
together$total_inserts = together$vj_insert + together$vd_insert + together$dj_insert

maf_dt = data.table()
for (snp in colnames(genotypes)[colnames(genotypes) != 'localID']){
    minor_allele = determine_true_minor_allele(snp, together)
    for (race in c(unique(together$race.g), 'Population')){
        maf = calculate_maf(snp, minor_allele, race, together)
        temp = data.table(snp = snp, ancestry_group = race, maf = maf, minor_allele = minor_allele)
        maf_dt = rbind(maf_dt, temp)
    }
}

bonferroni = 7.74e-9
sig_snps = find_significant_snps(c('vd_insert_no_pca', 'dj_insert_no_pca'), bonferroni)

sig_snps_together_mafs = maf_dt[snp %in% sig_snps$snp]


t_test = sig_snps_together_mafs %>%
    t_test(maf ~ ancestry_group, ref.group = 'Population')

ggboxplot(sig_snps_together_mafs, x = 'ancestry_group', y = 'maf', fill = 'ancestry_group', lwd = 1.5) +
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
    stat_pvalue_manual(t_test, label = 'p.adj.signif', size = 10, y.position = 12, remove.bracket = TRUE) +
    theme_classic() +
    theme(text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none")
    ggtitle('Minor Allele Frequency by Ancestry Group\nfor significant DNTT SNPs') +
    # geom_hline(data = average_all, aes(yintercept = population_maf_mean), size = 2.5, color = 'red', linetype = 2) +
    xlab('Ancestry Group') +
    ylab('Minor Allele Frequency') 

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/maf_by_race_boxplot.pdf'), plot = last_plot(), width = 15, height = 12, units = 'in', dpi = 750, device = 'pdf')


