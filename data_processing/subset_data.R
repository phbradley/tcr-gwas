library(data.table)
library(plyr)
library(readr)
library(stringr)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)

setwd("../_ignore/")

# ALL VGENE TRIMMING DATA--condensed by patient
condensed_trimming_data = as.data.table(read.table('emerson_stats/condensed_trim_data_all_patients.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE)[-1])

# subset data and write table for two patients
condensed_trimming_data_subset = condensed_trimming_data[patient_id == "HIP00110" | patient_id == "HIP00169"]
write.table(condensed_trimming_data_subset, file='subset_condensed_trimming_data_2_patients.tsv', quote=FALSE, sep='\t', col.names = NA)

# ALL SNP DATA
snps = snpgdsOpen("snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
snp_freq = snpgdsSNPRateFreq(snps, with.snp.id=TRUE)
snpgdsClose(snps)

# Filter snp data by missing rate less than 10%
snp_freq_dt = as.data.table(snp_freq)
snp_freq_dt_missing_rate_filtered = snp_freq_dt[MissingRate <= 0.1]

# write table for all filtered snp data
write.table(snp_freq_dt_missing_rate_filtered, file='snp_freq_data_missing_rate_filtered.tsv', quote=FALSE, sep='\t', col.names = NA)

# write table for subset of filtered snp data--100 snps
write.table(snp_freq_dt_missing_rate_filtered[1:100], file='subset_snp_freq_100_snps.tsv', quote=FALSE, sep='\t', col.names = NA)

