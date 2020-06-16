library(lme4)
library(data.table)

source("trimming_basic_regression_functions.R")


# Read in data

subset_snps = read.table("../_ignore/subset_snps_data_2500_snps.tsv", sep = "\t", fill=TRUE, header = TRUE)
trimming = as.data.table(read.table("../_ignore/condensed_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(trimming) = c("localID", "v_gene", "productive", "v_trim", "v_gene_count", "weighted_v_gene_count")

# Replace '3' genotypes with NA (missing genotype)
subset_snps[subset_snps==3]<-NA
subset_snps = as.data.table(subset_snps)

# Get ride of snp columns which have NA for all entries (missing for all individuals)
subset_snps_no_NA = subset_snps[,which(unlist(lapply(subset_snps, function(x)!all(is.na(x))))),with=F]

trimming_avg = trimming[, mean(v_trim), by = .(localID, productive)]
colnames(trimming_avg) = c("localID", "productive", "v_trim")

print("finished loading data")

# First, for each patient, we can create a regression to predict trimming length (averaging over all t_cells) given SNP minor allele frequency
# (here, averaging trimming length for all tcells for each patient)

#simple_regression_patient_average_results_productive = trimming_snp_regression_simple(subset_snps_no_NA, trimming_avg, productive = "True")

#print("finished simple_regression_patient_average_results_productive")

#simple_regression_patient_average_results_NOT_productive = trimming_snp_regression_simple(subset_snps_no_NA, trimming_avg, productive = "False")

#print("finished simple_regression_patient_average_results_NOT_productive")

# Next, for each patient, vgene we can create a regression to predict trimming length given SNP minor allele frequency
# (here, averaging trimming length for all tcells with the same vgene)

#simple_regression_patient_vgene_results_productive = trimming_snp_regression_simple(subset_snps_no_NA, trimming, productive = "True")

#print("finished simple_regression_patient_vgene_results_productive")

#simple_regression_patient_vgene_results_NOT_productive = trimming_snp_regression_simple(subset_snps_no_NA, trimming, productive = "False")

#print("finished simple_regression_patient_vgene_results_NOT_productive")


# Next, for each patient, vgene we can create a regression to predict trimming length given SNP minor allele frequency
# (here, averaging trimming length for all tcells with the same vgene) and weighting regression by proportion of t-cells with the specified vgenes, by patient

weighted_regression_patient_vgene_results_productive = trimming_snp_regression_weighted(subset_snps_no_NA, trimming, productive = "True")

print("finished weighted_regression_patient_vgene_results_productive")

weighted_regression_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(subset_snps_no_NA, trimming, productive = "False")

print("finished weighted_regression_patient_vgene_results_NOT_productive")

# Next, for each patient, vgene we can create a regression to predict trimming length given SNP minor allele frequency, but now, allowing for a subject variable intercept
# (here, averaging trimming length for all tcells with the same vgene)]

#simple_regression_patient_vgene_results_varying_int_subject = trimming_snp_regression_by_vgene_subject_varying_intercepts_subject(subset_snps_no_NA, trimming)

#print("finished simple_regression_patient_vgene_results_varying_int_subject")

## above does not work yet!!