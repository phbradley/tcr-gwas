
source("/home/mrussel2/tcr-gwas/trimming_regression/run_bootstrap_regression_all_snps_functions_cluster.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/make_snp_file.R")

library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

count = 10

all_snps = open_regressed_file_and_subset_by_pval(significance_cutoff = 5e-5, trim_type = args[2], random_effects = 'True', maf_cutoff = 0.05)

snp_subset = make_snp_file_subset_by_count_and_index(all_snps, count = count, index = as.numeric(args[1]))

if (nrow(snp_subset) != 0){
    genotype_subset = make_genotype_file_given_random_snps(snp_subset)

    print(paste0('finished compiling data for index ', args[1], " of ", nrow(all_snps)/count))

    # Run regression/bootstrap
    run_snps_trimming_snp_list_cluster(snp_list = snp_subset, genotype_list = genotype_subset, trim_type = args[2], gene_type = 'same', condensing = 'by_gene', gene_conditioning = 'True', weighting = 'True', random_effects = 'True', repetitions = as.numeric(args[4]), write_table = 'True', ncpus = as.numeric(args[3]), d_infer = 'True', maf_cutoff = 'False')
    print(paste0("Finished regressions for index ", get(args[1]), " of ", nrow(all_snps)/as.numeric(get(args[1]))))
} else {
    print('finished--no regressions necessary')
}
