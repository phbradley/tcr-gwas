set_regression_parameters <- function(phenotype){
    source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/', phenotype, '.R'))
    parameters = c(phenotype = phenotype, condensing_variable = CONDENSING_VARIABLE, keep_missing_d_gene = KEEP_MISSING_D_GENE, pca_count = as.character(PCA_COUNT))
    return(parameters)
}
