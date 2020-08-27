source("run_bootstrap_regression_all_snps_functions.R", chdir = TRUE)

# Read in snp list
snp_list = as.data.table(read.table("../_ignore/phil_sig_snps/top100_trimming_pvals_for_maggie.tsv", sep = "\t", fill=TRUE, header = TRUE))[-c(1:4)]
snps_list = unique(snp_list[chromosome == 10]$snpid)

# Run simple weighted and unweighted regression using Phil's condensed data (i.e. no d gene inference)

# read and reshape phil's trimming data
phil_trimming = as.data.table(read.table("../../_ignore/tmp.total_trims.tsv", sep = "\t", fill=TRUE, header = TRUE))
phil_IF = phil_trimming[,-c("v_trim_OF", 'j_trim_OF', 'd0_trim_OF', 'd1_trim_OF', 'total_trim_IF', 'total_trim_OF')]
colnames(phil_IF) = c("localID", 'v_trim', 'd0_trim', 'd1_trim', 'j_trim')
phil_IF$productive = as.logical('TRUE')
phil_OF = phil_trimming[,-c("v_trim_IF", 'j_trim_IF', 'd0_trim_IF', 'd1_trim_IF', 'total_trim_IF', 'total_trim_OF')]
colnames(phil_OF) = c("localID", 'v_trim', 'd0_trim', 'd1_trim', 'j_trim')
phil_OF$productive = as.logical('FALSE')
phil_trimming_reshaped = rbind(phil_OF, phil_IF)

assign('trimming_data', as.data.table(read.table(paste0("../_ignore/by_patient_condensed_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
setnames(trimming_data, "patient_id", "localID")

phil_trimming_reshaped = merge(phil_trimming_reshaped, trimming_data[,c('localID', 'productive', 'tcr_count')])
write.table(phil_trimming_reshaped, file= "../_ignore/by_patient_condensed_data_all_patients_phil.tsv", quote=FALSE, sep='\t', col.names = NA)

for (trim in c('d0_trim', 'd1_trim')){
    assign(paste0('by_gene_weighted_', trim), run_snps_trimming_snp_list(snp_id_list = snps_list, trim_type = trim, gene_type = 'd_gene', condensing = 'by_gene', gene_conditioning = 'True', weighting = 'True', repetitions = 100, write_table = 'False'))
    for (condense in c('by_patient', 'phil')){
        for (weight in c('True', 'False')){
            weight_name = ifelse('True', '_weighted_', '_no_weight_')
            assign(paste0(condense, weight_name, trim), run_snps_trimming_snp_list(snp_id_list = snps_list, trim_type = trim, gene_type = 'd_gene', condensing = condense, gene_conditioning = 'True', weighting = weight, repetitions = 100, write_table = 'False'))
        }
    }
}

together = data.table()
for (trim in c('d0_trim', 'd1_trim')){
    temp1 = data.table()
    temp1 = get(paste0('by_gene_weighted_', trim))
    temp1$condensing = 'by_gene'
    temp1$weighting = 'True'
    temp1$trim_type = trim
    together = rbind(together, temp1)
    for (condense in c('by_patient', 'phil')){
        for (weight in c('True', 'False')){
            temp2 = data.table()
            weight_name = ifelse('True', '_weighted_', '_no_weight_')
            temp2 = get(paste0(condense, weight_name, trim))[,-c('pval_py')]
            temp2$condensing = condense
            temp2$weighting = weight
            temp2$trim_type = trim
            together = rbind(together, temp2)
        }
    }
}

