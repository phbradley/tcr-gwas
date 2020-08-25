library("data.table")
library("plyr")
library("readr")
library("stringr")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")


compile_genotype_data <- function(snps_gds_file, snp_id){
    # get snp genotype
    genotype = snpgdsGetGeno(snps_gds, snp.id=snp, with.id = TRUE)
    
    # combine genotype, sample ids
    genotypes_dt = data.table(as.numeric(as.character(genotype$sample.id)), genotype$genotype)
    colnames(genotypes_dt) = c("scanID", paste0("snp",snp))
    # Convert subject names : 
    subject_id_mapping = as.data.table(read.table('../_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
    snps_genotypes = merge(genotypes_dt, subject_id_mapping, by = "scanID")
    return(snps_genotypes[,-c('scanID')])
}

filter_by_productivity <- function(condensed_trimming_dataframe, productive){
    # subset trimming data to include only productive or not productive entires
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "TRUE"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "FALSE"]
    } 
    return(condensed_trimming_dataframe)
}