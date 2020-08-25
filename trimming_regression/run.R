
source("run_bootstrap_regression_all_snps_functions.R")
source("find_correlated_snps.R")
source("make_snp_file.R")

args = commandArgs(trailingOnly=TRUE)


# Read in snp list
snp_list = as.data.table(read.table("../_ignore/phil_sig_snps/top100_trimming_pvals_for_maggie.tsv", sep = "\t", fill=TRUE, header = TRUE))[-c(1:4)]
artemis = snp_file(chromosome = 10, position1= 14397359, position2= 15454432)
tdt = snp_file(chromosome = 10, position1= 95804409, position2= 96838564)
mhc = snp_file(chromosome = 6, position1= 0, position2= 4050069)
rag = snp_file(chromosome = 11, position1= 36010709, position2= 37093156)

snp_data = get(args[2])
print(snp_data)

# Run regression/bootstrap
run_snps_trimming_snp_list(snp_id_list = unique(snp_data$snp), trim_type = args[1], condensing = 'by_patient', gene_conditioning = 'True', weighting = 'True', repetitions = 0, write_table = 'True')
print(paste0("Finished regressions for ", args[1], " for ", args[2])
#run_snps_trimming_snp_list(snp_id_list = unique(snp_list$snpid), trim_type = args[1], condensing = 'by_patient', gene_conditioning = 'False', weighting = 'False', repetitions = 100)

# Run simple regression/bootstrap
#run_snps_trimming_snp_list(snp_id_list = unique(snp_list$snpid), trim_type = args[1], varying_int = "False", repetitions = 100)


# Get correlated snp 
#chr7_data = as.data.table(read.table("../_ignore/tcrb_vgene_oof_hg19_pvals.tsv", header = TRUE))

# find correlated snps
#group_p_values(snp_id_list = unique(snp_list$snpid), snp_meta_data = snp_list, productivity = "productive", trim_type = args[1], bonferroni, proximity_cutoff, pvalue_cutoff)
#group_p_values(snp_id_list = unique(snp_list$snpid), snp_meta_data = snp_list, productivity = "NOT_productive", trim_type = args[1], bonferroni, proximity_cutoff, pvalue_cutoff)