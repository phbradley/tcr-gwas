library("tidyverse")
library("SNPRelate")


compile_genotype_data <- function(snps_gds_file, snp_id){
    # get snp genotype
    genotype = snpgdsGetGeno(snps_gds_file, snp.id=snp_id, with.id = TRUE)
    
    # combine genotype, sample ids
    genotypes_dt = data.frame(as.numeric(as.character(genotype$sample.id)), genotype$genotype)
    colnames(genotypes_dt) = c("scanID", paste0("snp",snp_id))
    # Convert subject names : 
    subject_id_mapping = as.data.frame(read.table('/home/mrussel2/tcr-gwas/_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
    snps_genotypes = merge(genotypes_dt, subject_id_mapping, by = "scanID")
    return(snps_genotypes[,c('localID', paste0("snp",snp_id))])
}

filter_by_productivity <- function(condensed_trimming_dataframe, productive){
    # subset trimming data to include only productive or not productive entires
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe %>% filter(productive == "TRUE")
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe %>% filter(productive == "FALSE")
    } 
    return(condensed_trimming_dataframe)
}

compile_trimming_data_cross <- function(){
    trimming_data_by_gene_all = data.frame()
    for (trim in c('v_trim', 'd0_trim', 'j_trim')){
        assign(paste0(trim, '_trimming_data'), as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/condensed_", trim, "_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        setnames(get(paste0(trim, '_trimming_data')), "patient_id", "localID")
        setnames(get(paste0(trim, '_trimming_data')), paste0(substring(trim, 1, 1), '_gene'), "gene")
        setnames(get(paste0(trim, '_trimming_data')), paste0(substring(trim, 1, 1), '_gene_count'), "gene_count")
        setnames(get(paste0(trim, '_trimming_data')), paste0('weighted_', substring(trim, 1, 1), '_gene_count'), "weighted_gene_count")
        gene_type = data.frame(gene_class = rep(paste0(substr(trim, 1, 1), '_gene'), nrow(get(paste0(trim, '_trimming_data')))))
        assign(paste0(trim, '_trimming_data'), cbind(get(paste0(trim, '_trimming_data')), gene_type))
        trimming_data_by_gene_all = rbind(trimming_data_by_gene_all, get(paste0(trim, '_trimming_data')))
    }
    return(trimming_data_by_gene_all)
}

# remove snp_genotype columns that are either all NA, or only have one genotype (for everyone)
remove_matrix_column_by_genotype <- function(genotype_matrix){
    for (snp in colnames(genotype_matrix)){
        genotypes = unique(genotype_matrix[,snp])
        nonNA_genotypes = genotypes[genotypes != 3]
        if (length(nonNA_genotypes) <= 1){
            genotype_matrix = genotype_matrix[, colnames(genotype_matrix) != snp]
        }
    }
    # replace 3 with NA
    genotype_matrix[genotype_matrix == 3] <- NA
    return(genotype_matrix)
}
