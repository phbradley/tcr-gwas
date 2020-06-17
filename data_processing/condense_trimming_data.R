library(data.table)
library(gdsfmt)
library(plyr)
library(readr)
library(stringr)

setwd("../_ignore/emerson_stats/")
dir = getwd()
files = list.files(path=dir, pattern="*.tsv", full.names=TRUE)

condense_trim_data <- function(gene_type){
    condensed_trim_data = data.table()
    for (file in files){
        temp_file = data.table()
        temp_file_condensed = data.table()
        temp_file = read.table(file, sep = "\t", fill=TRUE, header = TRUE)
        temp_file = as.data.table(temp_file)
        file_name = str_split(file, "/")[[1]][7]
        file_root_name = str_split(file_name, ".tsv")[[1]][1]
        patient_id = str_split(file_root_name, "_")[[1]][3]
        temp_file$patient_id = patient_id
        if (gene_type == 'v_gene'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, v_gene, productive)]
            temp_file_condensed = temp_file[,mean(v_trim), by = .(patient_id, v_gene, productive)]
            together = merge(temp_file_condensed,temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "v_gene", "productive", "v_trim", "v_gene_count")
            together$weighted_v_gene_count = together$v_gene_count/nrow(temp_file)
        } else if (gene_type == 'd_gene'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, d_gene, productive)]
            temp_file_condensed_d0 = temp_file[,mean(d0_trim), by = .(patient_id, d_gene, productive)]
            colnames(temp_file_condensed_d0) = c("patient_id", "d_gene", "productive", "d0_trim")
            temp_file_condensed_d1 = temp_file[,mean(d1_trim), by = .(patient_id, d_gene, productive)]
            colnames(temp_file_condensed_d1) = c("patient_id", "d_gene", "productive", "d1_trim")
            d_together = merge(temp_file_condensed_d0, temp_file_condensed_d1, by = c("patient_id", "d_gene", "productive"))
            together= merge(d_together, temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "d_gene", "productive", "d0_trim", "d1_trim", "d_gene_count")
            together$weighted_d_gene_count = together$d_gene_count/nrow(temp_file)
        } else if (gene_type == 'j_gene'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, j_gene, productive)]
            temp_file_condensed = temp_file[,mean(j_trim), by = .(patient_id, j_gene, productive)]
            together = merge(temp_file_condensed,temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "j_gene", "productive", "j_trim", "v_gene_count")
            together$weighted_j_gene_count = together$j_gene_count/nrow(temp_file)
        } 
        condensed_trim_data = rbind(condensed_trim_data, together)
        print(paste0(file, "processed for ", gene_type))
    }
    print(paste0("finished with ", gene_type))
    return(condensed_trim_data)
}

#v_gene_trimming = condense_trim_data('v_gene')
#write.table(v_gene_trimming, file='../condensed_v_trim_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)


d_gene_trimming = condense_trim_data('d_gene')
write.table(d_gene_trimming, file='../condensed_d_trim_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)

j_gene_trimming = condense_trim_data('j_gene')
write.table(j_gene_trimming, file='../condensed_j_trim_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)
