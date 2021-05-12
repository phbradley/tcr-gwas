library(data.table)
library(plyr)

source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/manhattan_plot_functions.R'))

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotations.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/analysis_scripts/analysis_functions.R'))


trimming_data = compile_manhattan_plot_data(c('v_trim', 'j_trim', 'd1_trim', 'd0_trim'))

artemis_snps = fread(paste0(PROJECT_PATH, '/tcr-gwas/_ignore/snp_data/artemis_missense_snps.tsv'))

artemis_hits = trimming_data[snp %in% artemis_snps$snpID]

print(artemis_snps[snpID %in% unique(artemis_hits$snp)])
print(artemis_hits)
