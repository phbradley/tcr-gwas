library(data.table)
setDTthreads(1)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(cowplot)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

features = c('gene_usage')

dataframe = compile_manhattan_plot_data(features) 
file_name = paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/gene_usage_top_association_distributions')

max_associations_by_gene = dataframe[, min(pvalue, na.rm = TRUE), by = .(productive, gene)]
colnames(max_associations_by_gene) = c('productive', 'gene', 'min_pvalue')
max_associations_by_gene[min_pvalue == Inf, min_pvalue := NA]

max_associations_by_gene[, gene_type := paste0(substring(gene, 4, 4), '-gene')]
max_associations_by_gene[productive == TRUE, productivity := 'productive']
max_associations_by_gene[productive == FALSE, productivity := 'non-productive']

significance_cutoff = determine_significance_cutoff(0.05, type = 'genome-wide', dataframe, phenotype = unique(dataframe$gene))

max_associations_by_gene = max_associations_by_gene[order(-gene_type)]
feature_colors = brewer.pal(n = max(4, length(unique(max_associations_by_gene$gene_type))), name = "Set2")
names(feature_colors) = unique(max_associations_by_gene$gene_type)

feature_colors = feature_colors[names(feature_colors) %in% unique(max_associations_by_gene$gene_type)]

max_associations_by_gene$gene_type = factor(max_associations_by_gene$gene_type, levels = c('V-gene', 'D-gene', 'J-gene'))
plot = ggplot(max_associations_by_gene, aes(x = -log10(min_pvalue))) +
    facet_grid(productivity ~ gene_type) +
    geom_histogram(aes(y = stat(width*density), group = gene_type, fill = gene_type), position = 'identity', binwidth = 3) +
    geom_vline(xintercept = -log10(significance_cutoff), color = 'grey70', size=2.5) +
    scale_fill_manual(values = feature_colors) 

final_plot = plot + theme_cowplot() + theme(legend.position = 'none', panel.spacing=unit(1.5, "lines"), text = element_text(size = 35), axis.text = element_text(size = 24), axis.line = element_blank(), axis.ticks = element_line(color = 'gray60', size = 1.5)) + coord_cartesian(clip="off") + ggtitle('') + background_grid(major = 'xy') + panel_border(color = 'gray60', size = 1.5) + ylab('Proportion of genes') + xlab('-log10(p-value)')  


ggsave(paste0(file_name, '.pdf'), plot = final_plot, width = 25, height = 18, units = 'in', dpi = 750, device = 'pdf')
saveRDS(final_plot, file = paste0(file_name, '.rds'))

