
source("run_bootstrap_regression_all_snps_functions.R")
source("correlate_snps.R")

args = commandArgs(trailingOnly=TRUE)

#run_snps_trimming(100, 5, 'v_trim')
#run_snps_trimming(100, 5, 'd0_trim')
#run_snps_trimming(100, 5, 'd1_trim')
#run_snps_trimming(100, 5, 'j_trim')

# Read in snp list
snp_list = as.data.table(read.table("../_ignore/tcrb_vgene_oof_hg19_pvals.tsv", sep = "\t", fill=TRUE, header = TRUE))
snp_list_temp = snp_list[500:505]

# Run regression/bootstrap for v_gene
run_snps_trimming_snp_list(snp_id_list = unique(snp_list_temp$snpnum), trim_type = args[1], varying_int = "TRUE", repetitions = 100)

# Get correlated snp 
chr7_data = as.data.table(read.table("../_ignore/tcrb_vgene_oof_hg19_pvals.tsv", header = TRUE))


correlated_snps_productive = group_p_values(snp_meta_data = chr7_data, productivity = "productive", trim_type = args[1], bonferroni, proximity_cutoff, pvalue_cutoff)
correlated_snps_NOT_productive = group_p_values(snp_meta_data = chr7_data, productivity = "NOT_productive", trim_type = args[1], bonferroni, proximity_cutoff, pvalue_cutoff)