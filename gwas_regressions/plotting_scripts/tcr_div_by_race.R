library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

tcr_div = fread(file = paste0(PROJECT_PATH, '/tcr-gwas/_ignore/emerson_tcrdiv_scores.tsv'))
setnames(tcr_div, 'hipid', 'localID')
setnames(tcr_div, 'tcrdiv', 'tcr_div')
ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))[,c('localID', 'race.g')]

together = merge(tcr_div, ethnicity, by = 'localID')[,c('localID', 'tcr_div', 'race.g')]
average_all = data.table(tcr_div_mean = mean(together$tcr_div))
ggplot(together, aes(x=race.g, y=tcr_div, fill = race.g)) +
    geom_boxplot(lwd = 1.5) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 4, alpha = 0.5) +
    stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", size = 10) +
    theme_classic() + 
    theme(text = element_text(size = 40), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
    ggtitle('TCR Repertoire Diversty by Race') +
    geom_hline(data = average_all, aes(yintercept = tcr_div_mean), size = 2.5, color = 'red', linetype = 2) +
    xlab('Racial Group') +
    ylab('TCR Diversity') 

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/tcr_div_by_race.pdf'), plot = last_plot(), width = 15, height = 12, units = 'in', dpi = 750, device = 'pdf')

