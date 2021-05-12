GENE_TYPE <<- 'd_gene'
PHENOTYPE_CLASS <<- 'allele_status_trimming'

CONDENSING_VARIABLE <<- 'by_gene_cdr3'
CONDITIONING_VARIABLE <<- 'cdr3_gene_group'
KEEP_MISSING_D_GENE <<- 'False'
REPETITIONS <<- 100
BOOTSTRAP_PVALUE_CUTOFF <<- 5e-5
PCA_COUNT <<- 8

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/phenotype_classes/', PHENOTYPE_CLASS, '_class_functions.R'))

