source('/home/mrussel2/tcr-gwas/plotting_scripts/insertion_distributions_by_gene.R')

plot_insertions_distribution_by_gene(c('v_gene', 'd_gene', 'j_gene'), c('vj_insert', 'vd_insert', 'dj_insert'), by_allele = 'True', frequency_filter = 0.01)
