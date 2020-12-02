# This script takes four arguments (1: starting snp position, 2: trim type, 3:
# number of cpus, 4: project_path, 5: output_path) 
args = commandArgs(trailingOnly=TRUE)


gene_name = args[1]
stopifnot(gene_name %in% c('dntt')
TRIM_TYPE= args[2]

# restrict threads
library(data.table)
setDTthreads(1)
library(ggplot2)
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

file_name =paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/by_race/', gene_name, '_', TRIM_TYPE, '_pc_by_asian_comparison.tsv') 

data = fread(file_name)

plot = ggplot(data, aes(x = -log10(pvalue_without_correction), y = -log10(pvalue_with_correction), color = pca_correction, shape = productivity)) +
    geom_point(size = 3, alpha = 0.5) +
    # geom_smooth(method = lm) +
    geom_abline(size = 2) + 
    ggtitle(paste0('P-value comparison at ', gene_name, ' locus for ', TRIM_TYPE)) +
    theme_classic() + 
    theme(text = element_text(size=16)) 


ggsave(paste0('figures/', gene_name, '_', TRIM_TYPE, '_pc_by_asian_comparison.png'), plot = plot, height=10, width=10, units="in", dpi = 500)
