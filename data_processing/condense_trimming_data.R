library(data.table)
library(gdsfmt)
library(plyr)
library(readr)
library(stringr)

setwd("../_ignore/emerson_stats/")
dir = getwd()
files = list.files(path=dir, pattern="*.tsv", full.names=TRUE)

condensed_trim_data = data.table()
for (file in files[2]){
    temp_file = data.table()
    temp_file_condensed = data.table()
    temp_file = read.table(file, sep = "\t", fill=TRUE, header = TRUE)
    temp_file = as.data.table(temp_file)
    file_name = str_split(file, "/")[[1]][7]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    patient_id = str_split(file_root_name, "_")[[1]][3]
    temp_file$patient_id = patient_id
    temp_file_condensed = temp_file[,mean(v_trim), by = .(patient_id, v_gene, productive)]
    colnames(temp_file_condensed) = c("patient_id", "v_gene", "productive", "v_trim")
    condensed_trim_data = rbind(condensed_trim_data, temp_file_condensed)
    print(paste(file))
}

write.table(condensed_trim_data, file='../condensed_trim_data_all_patients.tsv', quote=FALSE, sep='\t', col.names = NA)