source("run_bootstrap_regression_all_snps_functions.R")


# Read in snp list
snp_list = as.data.table(read.table("../_ignore/phil_sig_snps/top100_trimming_pvals_for_maggie.tsv", sep = "\t", fill=TRUE, header = TRUE))[-c(1:4)]

snp_data = snp_list[11]

# Run regression/bootstrap
run_snps_trimming_snp_list(snp_id_list = unique(snp_data$snp), trim_type = 'v_trim', condensing = 'by_gene', gene_conditioning = 'True', weighting = 'True', repetitions = 0, write_table = 'False')