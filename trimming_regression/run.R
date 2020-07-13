
source("run_bootstrap_regression_all_snps_functions.R")
source("correlate_snps.R")

args = commandArgs(trailingOnly=TRUE)


# Read in snp list
snp_list = as.data.table(read.table("../_ignore/phil_sig_snps/top100_trimming_pvals_for_maggie.tsv", sep = "\t", fill=TRUE, header = TRUE))[-c(1:4)]

# Run regression/bootstrap for v_gene
run_snps_trimming_snp_list(snp_id_list = unique(snp_list$snpid), trim_type = args[1], varying_int = "True", repetitions = 100)

# Run regression/bootstrap for v_gene
run_snps_trimming_snp_list(snp_id_list = unique(snp_list$snpid), trim_type = args[1], varying_int = "False", repetitions = 100)


# Get correlated snp 
#chr7_data = as.data.table(read.table("../_ignore/tcrb_vgene_oof_hg19_pvals.tsv", header = TRUE))

# find correlated snps
#group_p_values(snp_id_list = unique(snp_list$snpid), snp_meta_data = snp_list, productivity = "productive", trim_type = args[1], bonferroni, proximity_cutoff, pvalue_cutoff)
#group_p_values(snp_id_list = unique(snp_list$snpid), snp_meta_data = snp_list, productivity = "NOT_productive", trim_type = args[1], bonferroni, proximity_cutoff, pvalue_cutoff)