GENE_TYPE <<- 'NA'
PHENOTYPE_CLASS <<- 'trimming_naive'

CONDENSING_VARIABLE <<- 'by_subject'
CONDITIONING_VARIABLE <<- FALSE
KEEP_MISSING_D_GENE <<- 'True'
REPETITIONS <<- FALSE
BOOTSTRAP_PVALUE_CUTOFF <<- FALSE
PCA_COUNT <<- 8
 
source(paste0(PROJECT_PATH, '/tcr-gwas/scripts/phenotype_functions/phenotype_classes/', PHENOTYPE_CLASS, '_class_functions.R'))

