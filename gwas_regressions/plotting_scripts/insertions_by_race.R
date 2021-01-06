library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)

source('config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

insertions = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/insertions_by_patient.tsv')[,-1]
setnames(insertions, 'patient_id', 'localID')
ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))[,c('localID', 'race.g')]

together = merge(insertions, ethnicity, by = 'localID')[,c('localID', 'vj_insert', 'vd_insert', 'dj_insert', 'race.g', 'productive')]
together = together[, lapply(.SD, mean), by = .(localID, productive, race.g)]
together$total_inserts = together$vj_insert + together$vd_insert + together$dj_insert
together$productive = ifelse(together$productive == TRUE, 'in-frame', 'out-of-frame')

average_all = together[, mean(total_inserts), by = .(productive)]
setnames(average_all, 'V1', 'total_insert_mean')

ggplot(together, aes(x=race.g, y=total_inserts, fill = race.g)) +
    facet_grid(cols = vars(productive)) +
    geom_boxplot(lwd = 1.5) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 4, alpha = 0.5) +
    stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", size = 10) +
    theme_classic() + 
    theme(text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
    ggtitle('Mean Total Insertions by Race') +
    geom_hline(data = average_all, aes(yintercept = total_insert_mean), size = 2.5, color = 'red', linetype = 2) +
    xlab('Racial Group') +
    ylab('Mean total inserts') 

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/insertions_by_race.pdf'), plot = last_plot(), width = 25, height = 12, units = 'in', dpi = 750, device = 'pdf')

