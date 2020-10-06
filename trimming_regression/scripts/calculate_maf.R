source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/make_snp_file.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_regression_data_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/manha_visualization.R")

library('RhpcBLASctl')
library('data.table')
setDTthreads(1)
omp_set_num_threads(1)
blas_set_num_threads(1)

library('foreach')
library('doParallel')

# This file creates MAF file

count = 10000

calculate_maf_by_snp_file <- function(genotype_data_filtered, complete_dt){
    genotype_dt = as.data.table(genotype_data_filtered)
    gt_individual_count = c()
    snp_list = c()
    for (column in colnames(genotype_dt)){
        gt_individual_count = c(gt_individual_count, sum(genotype_dt[[column]] %in% c(0,1,2)))
        snp_list = c(snp_list, column)
    }
    gt_individuals = data.table(snp = paste0('snp', snp_list), individual_count = gt_individual_count)
    
    gt_sums = data.table(snp = paste0('snp', colnames(genotype_data_filtered)), gt_sums = colSums(genotype_data_filtered, na.rm = TRUE))

    gt_stats = merge(gt_individuals, gt_sums)
    
    gt_stats$allele_freq = gt_stats$gt_sums/(2*gt_stats$individual_count)
    
    gt_stats$maf = ifelse(gt_stats$allele_freq < 0.5, gt_stats$allele_freq, 1-gt_stats$allele_freq)

    return(gt_stats[,c('snp', 'maf')])
}

registerDoParallel(cores=20)
foreach(snp_start = seq(1, 35481497, count)) %dopar% {
        RcppParallel::setThreadOptions(1L) 
        snp_data = snp_file_by_snp_start(snp_start = snp_start, count)
        genotype_data = compile_all_genotypes(snp_start = snp_start, count)
        genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)
        print(paste0('starting calculation for snps starting at ', snp_start))
        results = calculate_maf_by_snp_file(genotype_data_filtered)
        write.table(results, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_results/maf_snps_starting_', snp_start, '.tsv'), quote=FALSE, sep='\t', col.names = NA)
        }
stopImplicitCluster()

compile_all_maf_data()