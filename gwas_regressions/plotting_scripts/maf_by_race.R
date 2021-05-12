library(cowplot)
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
    for (race in c(unique(together$race.g), 'all')){
        maf = calculate_maf(snp, minor_allele, race, together)
        temp = data.table(snp = snp, ancestry_group = race, maf = maf, minor_allele = minor_allele)
        maf_dt = rbind(maf_dt, temp)
    }
}

bonferroni = 7.74e-9
sig_snps = find_significant_snps(c('vd_insert_no_pca', 'dj_insert_no_pca'), bonferroni)

sig_snps_together_mafs = maf_dt[snp %in% sig_snps$snp]

# reshape data
all_race_mafs = sig_snps_together_mafs[ancestry_group== 'all']
by_race_mafs = sig_snps_together_mafs[ancestry_group != 'all']
sig_snps_together_reshaped = merge(by_race_mafs, all_race_mafs, by = c('snp', 'minor_allele'))[,c('snp', 'ancestry_group.x', 'maf.x', 'maf.y')]
colnames(sig_snps_together_reshaped) = c('snp', 'ancestry_group', 'maf_by_race', 'maf_all_races')

plot = ggplot(sig_snps_together_reshaped) +
    geom_point(aes(x = maf_all_races, y = maf_by_race, color = ancestry_group), alpha = 0.5, size = 7) +
    geom_abline(slope = 1, intercept = 0,  size = 2) +
    # ggtitle('MAF by ancestry group versus by population \nfor DNTT significant snps') +
    theme_classic() +
    theme(text = element_text(size = 40)) +
    xlab('MAF computed across all individuals') +
    ylab('MAF computed by ancestry group') +
    labs(fill = "Ancestry Group") +
    scale_color_brewer(palette = 'Set2')

final_plot = plot + theme_cowplot() + theme(text = element_text(size = 40), axis.text = element_text(size = 24), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'xy') +  labs(color = 'Ancestry group')

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/maf_by_race.pdf'), plot = final_plot, width = 14, height = 10, units = 'in', dpi = 750, device = 'pdf')
