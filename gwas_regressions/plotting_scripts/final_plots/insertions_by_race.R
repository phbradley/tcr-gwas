library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(Cairo)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

insertions = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/insertions_by_patient.tsv')[,-1]
setnames(insertions, 'patient_id', 'localID')
ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))[,c('localID', 'race.g')]

together = merge(insertions, ethnicity, by = 'localID')[,c('localID', 'vj_insert', 'vd_insert', 'dj_insert', 'race.g', 'productive')]
together = together[, lapply(.SD, mean), by = .(localID, productive, race.g)]
together$total_inserts = together$vj_insert + together$vd_insert + together$dj_insert
together$productive = ifelse(together$productive == TRUE, 'productive', 'non-productive')

average_all = together[, mean(total_inserts), by = .(productive)]
setnames(average_all, 'V1', 'total_insert_mean')
setnames(together, 'race.g', 'ancestry_group')

t_test = together %>%
    group_by(productive) %>%
    t_test(total_inserts ~ ancestry_group, ref.group = 'all')

plot = ggboxplot(together, x = 'ancestry_group', y = 'total_inserts', fill = 'ancestry_group', lwd = 1.5) +
    facet_grid(cols = vars(productive)) +
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, alpha = 0.5) +
    stat_pvalue_manual(t_test, label = "p = {p}", size = 8, y.position = 12, remove.bracket = TRUE, family = 'Arial') +
    theme_classic() + 
    theme(text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
    ggtitle('Mean Total Insertions by Ancestry Group') +
    geom_hline(data = average_all, aes(yintercept = total_insert_mean), size = 2.5, color = 'red', linetype = 2) +
    xlab('Ancestry Group') +
    ylab('Mean total inserts') + 
    scale_fill_brewer(palette = 'Set2')

final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5, size = 24), axis.text.y = element_text(size = 24), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'xy') + panel_border(color = 'gray60', size = 1.5)
ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/insertions_by_race.pdf'), plot = final_plot, width = 25, height = 12, units = 'in', dpi = 750, device = cairo_pdf)

