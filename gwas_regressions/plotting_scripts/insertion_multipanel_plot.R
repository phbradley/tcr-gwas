library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

insertion_manhattan = readRDS(paste0(make_file_name('insertion'), '.rds'))
insertion_manhattan_dntt = readRDS(paste0(make_file_name('insertion', 'dntt'), '.rds'))

revised_dntt = insertion_manhattan_dntt + theme(legend.position = "none") + ggtitle('') + xlab('') + ylab('')
revised_insertion = insertion_manhattan + theme(legend.position = "none") + ggtitle('')
aligned_plots = align_plots(revised_insertion, revised_dntt, align = 'v', axis = 'l')
legend = get_legend(insertion_manhattan_dntt)

bottom_row = plot_grid(aligned_plots[[2]], legend, rel_widths = c(2, 1), nrow = 1, labels = c('B'), label_size = 40, align="none")

together = plot_grid(aligned_plots[[1]], 
                     bottom_row, 
                     ncol = 1, 
                     rel_heights = c(2, 1.2), 
                     label_size = 40,
                     labels = c('A'))

file_name = paste0(paste0(make_file_name('insertion'), '_with_dntt_multipanel.pdf'))
ggsave(file_name, plot = together, width = 25, height = 18, units = 'in', dpi = 750, device = 'pdf')

