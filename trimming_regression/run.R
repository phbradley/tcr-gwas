
source("run_bootstrap_regression_all_snps_functions.R")
args = commandArgs(trailingOnly=TRUE)

#run_snps_trimming(100, 5, 'v_gene')
#run_snps_trimming(100, 5, 'd0_gene')
#run_snps_trimming(100, 5, 'd1_gene')
#run_snps_trimming(100, 5, 'j_gene')

# Read in snp list
snp_list = as.data.table(read.table("../_ignore/tcrb_vgene_oof_hg19_pvals.tsv", sep = "\t", fill=TRUE, header = TRUE))

# Run regression/bootstrap for v_gene
run_snps_trimming_snp_list(snp_id_list = unique(snp_list$snpnum), gene_type = args[1])