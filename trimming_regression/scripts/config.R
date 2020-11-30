
#set global regression parameters based on trimming/insertion type
set_regression_parameters <- function(trim_type){
    type <- strsplit(trim_type, '_')[[1]][2]

    # Set regression parameters
    WEIGHTING <<- ifelse(type == 'insert', 'False', 'True')
    # also could condense 'by_gene', but this has more observations and, thus, covariates in the model
    CONDENSING <<- ifelse(type == 'insert', 'by_patient', 'by_cdr3')
    GENE_CONDITIONING <<- ifelse(type == 'insert', 'False', 'True')
    RANDOM_EFFECTS <<- 'False'
    D_INFER <<- ifelse(type == 'insert', 'False', 'True')
    REPETITIONS <<- ifelse(type == 'insert', 0, 100)
    BOOT_CUTOFF <<- 5e-5
    MAF_CUTOFF <<- 0.05
    PVALUE_VARIABLE <<- 'snp'
}

set_regression_parameters(TRIM_TYPE)
