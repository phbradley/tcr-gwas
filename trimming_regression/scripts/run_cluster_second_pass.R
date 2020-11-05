
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/make_snp_file.R")

library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

count = 10

# parse type (either 'insert' or 'trim')
stopifnot(args[2] %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'vj_insert'))
type = strsplit(args[2], '_')[[1]][2]


# set variables conditional on type (i.e. if we are doing insertion
# regressions, we want to condense the trimming data by patient without
# inferring d_gene, and regress not conditioning out TCRB gene and without 
# random effects)
condensing = ifelse(type == 'insert', 'by_patient', 'by_gene')
gene_conditioning = ifelse(type == 'insert', 'False', 'True')
random_effects = ifelse(type == 'insert', 'False', 'True')
d_infer = ifelse(type == 'insert', 'False', 'True')
pca_structure_correction = ifelse(args[5] == 'PCA', 'True', 'False')

condensing = ifelse(type == 'insert', 'by_patient', 'by_patient')
gene_conditioning = ifelse(type == 'insert', 'False', 'False')
random_effects = ifelse(type == 'insert', 'False', 'False')
d_infer = ifelse(type == 'insert', 'False', 'False')

all_snps = open_regressed_file_and_subset_by_pval(significance_cutoff = 5e-5, trim_type = args[2], random_effects, condensing, d_infer, maf_cutoff = 0.05)

snp_subset = make_snp_file_subset_by_count_and_index(all_snps, count = count, index = as.numeric(args[1]))


if (nrow(snp_subset) != 0){
    genotype_subset = make_genotype_file_given_random_snps(snp_subset)

    print(paste0('finished compiling data for index ', args[1], " of ", nrow(all_snps)/count))

    # Run regression/bootstrap
    run_snps_trimming_snp_list_cluster(snp_list = snp_subset, genotype_list = genotype_subset, trim_type = args[2], gene_type = 'same', condensing, gene_conditioning, weighting = 'True', random_effects, repetitions = as.numeric(args[4]), pca_structure_correction, write_table = 'True', ncpus = as.numeric(args[3]), d_infer, maf_cutoff = 'False')
    print(paste0("Finished regressions for index ", args[1], " of ", nrow(all_snps)/as.numeric(args[1])))
} else {
    print('finished--no regressions necessary')
}
