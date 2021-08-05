library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)
library(Cairo)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

manhattan_plot_gene(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
                    phenotype_values = c('V-D-gene junction N-insertions', 'D-J-gene junction N-insertions'),
                    plot_title = set_manhattan_plot_title('insertion', 'dntt'), 
                    file_name = make_file_name('insertion', 'dntt'), 
                    gene_subset = 'dntt')

manhattan_plot(dataframe = compile_manhattan_plot_data(c('vd_insert', 'dj_insert')), 
               phenotype_values = c('V-D-gene junction N-insertions', 'D-J-gene junction N-insertions'),
               plot_title = set_manhattan_plot_title('insertion'), 
               file_name = make_file_name('insertion'), 
               plotting_p_value_cutoff =  -log10(5e-3), 
               subsample_point_cutoff = -log10(0.001)) 

insertion_manhattan = readRDS(paste0(make_file_name('insertion'), '.rds')) + labs(color = NULL)
insertion_manhattan_dntt = readRDS(paste0(make_file_name('insertion', 'dntt'), '.rds')) + labs(color = NULL)


revised_insertion = insertion_manhattan + theme_cowplot(font_family = 'Arial') + theme(legend.direction = "horizontal",legend.position="bottom", legend.justification = 'left', panel.spacing.y=unit(3, "lines"), axis.ticks.x = element_blank(), panel.spacing.x=unit(0, "lines"), strip.placement = 'outside', strip.background.x = element_blank(), text = element_text(size = 40), axis.text.x = element_blank(), strip.text.x = element_text(size = 12), axis.text.y = element_text(size = 24), axis.line.y = element_line(color = 'gray60', size = 1.5), axis.line.x = element_blank(), axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + geom_hline(aes(yintercept = Inf), size = 1.5, color = 'gray60')+geom_hline(aes(yintercept = -Inf), size = 1.5, color = 'gray60') + scale_x_continuous() + geom_point(aes(x = Inf, y = -Inf), size = 3, color = 'gray60', shape = '\uFF5C') + ggtitle('') + ylim(2, 12.5) + background_grid(major = 'y')


revised_dntt = insertion_manhattan_dntt + theme_cowplot(font_family = 'Arial') + theme(legend.position = "none", panel.spacing.y=unit(3, "lines"), panel.spacing.x=unit(0, "lines"), strip.placement = 'outside', strip.background.x = element_blank(), text = element_text(size = 40), strip.text.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 24), axis.line.y = element_line(color = 'gray60', size = 1.5), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + geom_hline(aes(yintercept = Inf), size = 1.5, color = 'gray60')+geom_hline(aes(yintercept = -Inf), size = 1.5, color = 'gray60') + scale_x_continuous()  + ggtitle('') + xlab('Chr 10 Position\n(bp 97,864,085-97,898,321)') + ylim(2, 12.5) + background_grid(major = 'y')
   
theme_set(theme_cowplot(font_family = "Arial"))

aligned_plots = align_plots(revised_insertion, revised_dntt, align = 'h', axis = 'bt')

together = plot_grid(aligned_plots[[1]], aligned_plots[[2]], ncol = 2, rel_widths = c(2, 0.8), label_size = 40, labels = "AUTO", scale = 0.97)

file_name = paste0(paste0(make_file_name('insertion'), '_with_dntt_multipanel.pdf'))

ggsave(file_name, plot = together, width = 26, height = 18, units = 'in', dpi = 750, device = cairo_pdf)


