
source("run_bootstrap_regression_all_snps_functions_cluster.R")
source("find_correlated_snps.R")
source("make_snp_file.R")

args = commandArgs(trailingOnly=TRUE)


# Read in snp list
snp_data = snp_file_by_snp_start(snp_start = as.numeric(args[1]), 50000)

# Run regression/bootstrap
run_snps_trimming_snp_list_cluster(snp_id_list = unique(snp_data$snp), trim_type = args[2], gene_type = 'same', condensing = 'by_gene', gene_conditioning = 'True', weighting = 'True', repetitions = 100, write_table = 'True', numCores = 4)
print(paste0("Finished regressions for ", args[2], " for snps ", args[1], '-', as.character(as.numeric(args[1])+50000)))