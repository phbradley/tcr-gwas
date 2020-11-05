library("data.table")
library("gdsfmt")
library("plyr")
library("readr")
library("stringr")

args = commandArgs(trailingOnly=TRUE)

setwd("/home/mrussel2/tcr-gwas/_ignore/emerson_stats/")
dir = getwd()
files = list.files(path=dir, pattern="*.tsv", full.names=TRUE)

compile_trim_data <- function(){
    trim_by_gene_type_patient = data.table()
    count = 0
    # for each patient...
    for (file in files){
        count = count + 1
        # read file...
        temp_file = fread(file, sep = "\t", fill=TRUE, header = TRUE)
        file_name = str_split(file, "/")[[1]][7]
        file_root_name = str_split(file_name, ".tsv")[[1]][1]
        patient_id = str_split(file_root_name, "_")[[1]][3]
        temp_file$patient_id = patient_id
        trim_by_gene_type_patient = rbind(trim_by_gene_type_patient, data.table(patient_id = temp_file$patient_id, v_gene = temp_file$v_gene, d_gene = temp_file$d_gene, j_gene = temp_file$j_gene, v_trim = temp_file$v_trim, d0_trim = temp_file$d0_trim, d1_trim = temp_file$d1_trim, j_trim = temp_file$j_trim, productivity = temp_file$productive))
        print(paste0('finished ', count, ' of ', length(files)))
    }
    colnames(trim_by_gene_type_patient) = c('patient_id', 'v_gene', 'd_gene', 'j_gene', 'v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'productive')
    write.table(trim_by_gene_type_patient, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/_by_patient.tsv'), quote=FALSE, sep='\t', col.names = NA)
}

compile_trim_data()
