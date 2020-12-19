trimming = c(condensing_variable = 'by_gene_cdr3', infer_missing_d_gene = 'True', pca_count = '8')
trimming_naive = c(condensing_variable = 'by_subject', infer_missing_d_gene = 'True', pca_count = '8')

insertion = c(condensing_variable = 'by_subject', infer_missing_d_gene = 'False', pca_count = '8')
insertion_no_pca = c(condensing_variable = 'by_subject', infer_missing_d_gene = 'False', pca_count = '0')

p_addition_count = c(condensing_variable = 'by_gene_cdr3', infer_missing_d_gene = 'False', pca_count = '8')
p_addition_fraction = c(condensing_variable = 'by_gene_cdr3', infer_missing_d_gene = 'False', pca_count = '8') 

set_regression_parameters <- function(phenotype){
    if (phenotype %in% c('v_trim', 'j_trim', 'd0_trim', 'd1_trim')){
        parameters = c(phenotype = phenotype, trimming)
    } else if (phenotype %in% c('vd_insert', 'vj_insert', 'dj_insert')){
        parameters = c(phenotype = phenotype, insertion)
    } else if (phenotype %in% c('v_trim_naive', 'j_trim_naive', 'd0_trim_naive', 'd1_trim_naive')){
        parameters = c(phenotype = phenotype, trimming_naive)
    } else if (phenotype %in% c('vd_insert_no_pca', 'vj_insert_no_pca', 'dj_insert_no_pca')){
        parameters = c(phenotype = phenotype, insertion_no_pca)
    } else if (phenotype %in% c('v_pnucs_count', 'j_pnucs_count', 'd0_pnucs_count', 'd1_pnucs_count')){
        parameters = c(phenotype = phenotype, p_addition_count)
    } else if (phenotype %in% c('v_pnucs_fraction', 'j_pnucs_fraction', 'd0_pnucs_fraction', 'd1_pnucs_fraction')){ 
        parameters = c(phenotype = phenotype, p_addition_fraction)
    }
    return(parameters)
}
