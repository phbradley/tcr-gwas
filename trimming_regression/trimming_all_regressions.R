library(lme4)
library(data.table)

source("trimming_basic_regression_functions.R")
source("../data_processing/compile_snp_data_all.R")

print("finished processing snp data")

# Read in data
snps = genotypes_dt_bothID

trimming = as.data.table(read.table("../_ignore/condensed_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(trimming) = c("localID", "v_gene", "productive", "v_trim", "v_gene_count", "weighted_v_gene_count")

# Replace '3' genotypes with NA (missing genotype)
snps[snps==3]<-NA
snps = as.data.table(snps)

# Get ride of snp columns which have NA for all entries (missing for all individuals)
snps_no_NA = snps[,which(unlist(lapply(snps, function(x)!all(is.na(x))))),with=F]

print("finished loading data")


weighted_regression_patient_vgene_results_productive = trimming_snp_regression_weighted(snps_no_NA, trimming, productive = "True")
print("finished weighted_regression_patient_vgene_results_productive")

weighted_regression_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(snps_no_NA, trimming, productive = "False")
print("finished weighted_regression_patient_vgene_results_NOT_productive")

