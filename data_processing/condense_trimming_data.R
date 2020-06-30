library(data.table)
library(gdsfmt)
library(plyr)
library(readr)
library(stringr)

setwd("../_ignore/emerson_stats/")
dir = getwd()
files = list.files(path=dir, pattern="*.tsv", full.names=TRUE)
infer_d_gene <- function(patient_trimming_table){
    patient_trimming_table_NO_missing_d = patient_trimming_table[d_gene != '-']
    patient_trimming_table_missing_d = patient_trimming_table[d_gene == '-']
    triplet_counts = patient_trimming_table_NO_missing_d[,.N, by = .(v_gene, d_gene, j_gene)]
    for (row in 1:nrow(patient_trimming_table_missing_d)){
        v_gene_observed = patient_trimming_table_missing_d[row]$v_gene
        j_gene_observed = patient_trimming_table_missing_d[row]$j_gene

        empirical_dist = triplet_counts[v_gene == v_gene_observed & j_gene == j_gene_observed]
        if (nrow(empirical_dist) != 0) {
            empirical_dist$probability = empirical_dist$N/sum(empirical_dist$N)
            sample = rmultinom(1,1,empirical_dist$probability)
            index = which(sample == 1)
            patient_trimming_table_missing_d[row]$d_gene = empirical_dist[index]$d_gene
        } else {
            sample = rmultinom(1,1,rep(0.33333, 3))
            index = which(sample == 1)
            d_gene_options = unique(patient_trimming_table_NO_missing_d$d_gene)
            patient_trimming_table_missing_d[row]$d_gene = d_gene_options[index]
        } 
        if (patient_trimming_table_missing_d[row]$d_gene == "TRBD2*01" | patient_trimming_table_missing_d[row]$d_gene == "TRBD2*02"){
            patient_trimming_table_missing_d[row]$d0_trim = 4
            patient_trimming_table_missing_d[row]$d1_trim = 4
            patient_trimming_table_missing_d[row]$dj_insert = 4
            patient_trimming_table_missing_d[row]$vd_insert = 4
        } else if (patient_trimming_table_missing_d[row]$d_gene == "TRBD1*01"){
            patient_trimming_table_missing_d[row]$d0_trim = 3
            patient_trimming_table_missing_d[row]$d1_trim = 3
            patient_trimming_table_missing_d[row]$dj_insert = 3
            patient_trimming_table_missing_d[row]$vd_insert = 3
        }
    }
    inferred = rbind(patient_trimming_table_missing_d, patient_trimming_table_NO_missing_d)
    return(inferred)
}

condense_trim_data <- function(trim_type){
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
        temp_file = infer_d_gene(temp_file)
        if (trim_type == 'v_trim'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, v_gene, productive)]
            temp_file_condensed = temp_file[,mean(v_trim), by = .(patient_id, v_gene, productive)]
            together = merge(temp_file_condensed,temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "v_gene", "productive", "v_trim", "v_gene_count")
            together$weighted_v_gene_count = together$v_gene_count/nrow(temp_file)
        } else if (trim_type == 'd_trim'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, d_gene, productive)]
            temp_file_condensed_d0 = temp_file[,mean(d0_trim), by = .(patient_id, d_gene, productive)]
            colnames(temp_file_condensed_d0) = c("patient_id", "d_gene", "productive", "d0_trim")
            temp_file_condensed_d1 = temp_file[,mean(d1_trim), by = .(patient_id, d_gene, productive)]
            colnames(temp_file_condensed_d1) = c("patient_id", "d_gene", "productive", "d1_trim")
            d_together = merge(temp_file_condensed_d0, temp_file_condensed_d1, by = c("patient_id", "d_gene", "productive"))
            together= merge(d_together, temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "d_gene", "productive", "d0_trim", "d1_trim", "d_gene_count")
            together$weighted_d_gene_count = together$d_gene_count/nrow(temp_file)
        } else if (trim_type == 'j_trim'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, j_gene, productive)]
            temp_file_condensed = temp_file[,mean(j_trim), by = .(patient_id, j_gene, productive)]
            together = merge(temp_file_condensed,temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "j_gene", "productive", "j_trim", "j_gene_count")
            together$weighted_j_gene_count = together$j_gene_count/nrow(temp_file)
        } else if (trim_type == 'vj_insert'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, v_gene, j_gene, productive)]
            temp_file_condensed = temp_file[,mean(vj_insert), by = .(patient_id, v_gene, j_gene, productive)]
            together = merge(temp_file_condensed,temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "v_gene", "j_gene", "productive", "vj_insert", "vj_gene_count")
            together$weighted_vj_gene_count = together$vj_gene_count/nrow(temp_file)
        } else if (trim_type == 'dj_insert'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, d_gene, j_gene, productive)]
            temp_file_condensed = temp_file[,mean(dj_insert), by = .(patient_id, d_gene, j_gene, productive)]
            together= merge(temp_file_condensed, temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "d_gene", "j_gene", "productive", "dj_insert", "dj_gene_count")
            together$weighted_dj_gene_count = together$dj_gene_count/nrow(temp_file)
        } else if (trim_type == 'vd_insert'){
            temp_file_vgene_type_counts = temp_file[,.N, by = .(patient_id, v_gene, d_gene, productive)]
            temp_file_condensed = temp_file[,mean(vd_insert), by = .(patient_id, v_gene, d_gene, productive)]
            together = merge(temp_file_condensed,temp_file_vgene_type_counts)
            colnames(together) = c("patient_id", "v_gene", "d_gene", "productive", "vd_insert", "vd_gene_count")
            together$weighted_vd_gene_count = together$vd_gene_count/nrow(temp_file)
        } 
        condensed_trim_data = rbind(condensed_trim_data, together)
        print(paste0(file, "processed for ", trim_type))
    }
    print(paste0("finished with ", trim_type))
    return(condensed_trim_data)
}

#v_gene_trimming = condense_trim_data('v_trim')
#write.table(v_gene_trimming, file='../condensed_v_trim_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)


#d_gene_trimming = condense_trim_data('d_trim')
#write.table(d_gene_trimming, file='../condensed_d_trim_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)

#j_gene_trimming = condense_trim_data('j_trim')
#write.table(j_gene_trimming, file='../condensed_j_trim_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)

vj_gene_insert = condense_trim_data('vj_insert')
write.table(vj_gene_insert, file='../condensed_vj_insert_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)

dj_gene_insert = condense_trim_data('dj_insert')
write.table(dj_gene_insert, file='../condensed_dj_insert_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)

vd_gene_insert = condense_trim_data('vd_insert')
write.table(vd_gene_insert, file='../condensed_vd_insert_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)