library(GWASTools)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(GGally)
library(Cairo)

source('config/config.R')
source('config/file_paths.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))[,c('localID', 'scanID', 'race.g')]

counts = ethnicity[, .N, by = race.g]
counts$race.g = factor(counts$race.g, levels = c('Caucasian', 'Hispanic', 'African', 'Asian', 'Native American', 'Middle Eastern'))

plot = ggplot(counts) +
    geom_bar(aes(x = race.g, y = N, fill = race.g), stat = 'identity') +
    xlab('Ancestry group') +
    ylab('Proportion of sample population') +
    scale_fill_brewer(palette = 'Set2')

final_plot = plot + theme_cowplot(font_family = 'Arial') + theme(legend.position = 'none', text = element_text(size = 25), axis.text = element_text(size = 15), axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15), axis.line = element_blank(), axis.ticks = element_blank()) + background_grid(major = 'y') 

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/population_composition.pdf'), plot = final_plot, width = 8, height = 8, units = 'in', dpi = 750, device = cairo_pdf)


