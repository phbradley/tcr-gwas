source('config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/file_paths.R'))

library(data.table)
setDTthreads(1)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

phenotype_list <<- args[1]

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/src/regression_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plot_src/plotting_functions/manhattan_plot_functions.R'))

dataframe = compile_manhattan_plot_data(phenotype_list)

bonferroni = 0.05/length(unique(dataframe$snp))
print(paste0('The bonferroni cutoff is ', bonferroni))

sigs = dataframe[pvalue < bonferroni]
sigs = sigs[order(pvalue)]
print(paste0('There are ', nrow(sigs), ' total significant associations'))

genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
chr = c(10, 6, 10, 11, 7, 14)
pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

gene_annotations = data.table(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)

artemis = gene_annotations[genes == 'artemis']


artemis_sigs = sigs[chr == artemis$chr & hg19_pos < artemis$pos2 & hg19_pos > artemis$pos1]

print(paste0('There are ', nrow(artemis_sigs), ' total significant associations within the artemis locus')) 

tcrb = gene_annotations[genes == 'tcrb']

tcrb_sigs = sigs[chr == tcrb$chr & hg19_pos < tcrb$pos2 & hg19_pos > tcrb$pos1]

print(paste0('There are ', nrow(tcrb_sigs), ' total significant associations within the tcrb locus')) 

dntt = gene_annotations[genes == 'dntt']

dntt_sigs = sigs[chr == dntt$chr & hg19_pos < dntt$pos2 & hg19_pos > dntt$pos1]

