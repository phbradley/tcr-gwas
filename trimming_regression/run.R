
source("run_bootstrap_regression_all_snps_functions.R")

run_snps_trimming(100, 5, 'v_gene')
run_snps_trimming(100, 5, 'd0_gene')
run_snps_trimming(100, 5, 'd1_gene')
run_snps_trimming(100, 5, 'j_gene')

#run_snps_trimming_parallel(1, 1000, 'v_gene')