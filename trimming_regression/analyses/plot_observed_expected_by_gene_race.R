# This script takes four arguments (1: starting snp position, 2: trim type, 3:
# number of cpus, 4: project_path, 5: output_path) 
args = commandArgs(trailingOnly=TRUE)


gene_name = args[1]
stopifnot(gene_name %in% c('dntt'))

trimming_type = args[2]
pca = args[3]
ncpu = args[4]
project_path = args[5]
output_path = args[6]
pcatype = args[7]

# restrict threads
library(data.table)
setDTthreads(1)
library(ggplot2)
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/by_race/', gene_name, '_', trimming_type, '_pca-', pca, '_', pcatype, '_by_race.tsv')

pvals = fread(file = file_name)[,-c(1)]

pvals_all_races = pvals[race == 'all_together']
pvals_by_race = pvals[race != 'all_together' & pvalue != 'NaN' & pvalue != 'NA']

pvals_all_races_small = pvals_all_races[,c('snp', 'pvalue', 'productivity')]
colnames(pvals_all_races_small) = c('snp', 'all_race_pvalue', 'productivity')

together = merge(pvals_by_race, pvals_all_races_small, by = c('snp', 'productivity'))

plot = ggplot(together) +
    geom_point(aes(x = -log10(all_race_pvalue), y = -log10(pvalue), color = race, shape = productivity), size = 3, alpha = 0.5) +
    geom_abline(size = 2) + 
    ggtitle(paste0('P-value comparison at ', gene_name, ' locus for ', trimming_type)) +
    theme_bw() + 
    theme(text = element_text(size=16)) +
    xlab('Population structure corrected -log10(pvalues) for all races combined') +
    ylab('Un-corrected -log10(pvalues) by race')


ggsave(paste0('figures/', gene_name, '_', trimming_type, '_pca-', pca, '_', pcatype, '_by_race.png'), plot = plot, height=10, width=10, units="in", dpi = 500)
