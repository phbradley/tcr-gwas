GENE_TYPE <<- 'NA' 
CONDENSING_VARIABLE <<- 'by_subject'
CONDITIONING_VARIABLE <<- FALSE
KEEP_MISSING_D_GENE <<- 'False'
REPETITIONS <<- FALSE
BOOTSTRAP_PVALUE_CUTOFF <<- FALSE
PCA_COUNT <<- 8
 
PHENOTYPE_CLASS <<- 'missing_d_fraction_j1_subset'

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/phenotype_classes/', PHENOTYPE_CLASS, '_class_functions.R'))

