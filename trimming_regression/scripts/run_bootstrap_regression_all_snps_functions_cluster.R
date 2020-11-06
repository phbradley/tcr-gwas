library('plyr')
library("tidyverse")
library('foreach')
library('doParallel')
library('RhpcBLASctl')
library('stringr')
library('data.table')
setDTthreads(threads = 1)
omp_set_num_threads(1)
blas_set_num_threads(1)


source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/regression_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/bootstrap_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_data_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/execute_regression_function.R")

# This function runs regressions on the cluster
run_snps_trimming_snp_list_cluster <- function(snp_list, genotype_list, trim_type, pca_structure_correction, pvalue_boot_threshold, write_table, ncpus, maf_cutoff, data_file_path){
 
    # parse type (either 'insert' or 'trim')
    stopifnot(args[2] %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'vj_insert'))
    type = strsplit(args[2], '_')[[1]][2]
 
    # Set regression parameters
    weighting = 'True'
    condensing = ifelse(type == 'insert', 'by_patient', 'by_gene')
    gene_conditioning = ifelse(type == 'insert', 'False', 'True')
    random_effects = ifelse(type == 'insert', 'False', 'True')
    d_infer = ifelse(type == 'insert', 'False', 'True')
    repetitions = ifelse(type == 'insert', 0, 100)
    
    regression_dataframe = data.frame()

    # import condensed trimming file
    trimming_data = compile_condensed_trimming_data(trim_type, d_infer, condensing)

    #filter out snps with maf below cutoff
    if (maf_cutoff == "False"){
        list_of_snps = as.numeric(colnames(genotype_list)[colnames(genotype_list) != 'localID'])
    } else {
        maf_file_name = '/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'
        if (!file.exists(maf_file_name)){
            system(command = paste0("Rscript /home/mrussel2/tcr-gwas/trimming_regression/scripts/calculate_maf.R ", ncpus))
        }
        maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]
        maf_data_filtered = maf_data %>% filter(maf >= maf_cutoff)
        maf_data_snps = gsub('snp', '', maf_data_filtered$snp)
        list_of_snps = as.numeric(intersect(colnames(genotype_list), maf_data_snps))
    }

    count = 0
    registerDoParallel(cores=ncpus)
    results = foreach(snp=list_of_snps, .combine='rbind') %dopar% {
        RcppParallel::setThreadOptions(1L) 
        execute_regression(snp, list_of_snps, snp_list, genotype_list, trim_type, condensing, trimming_data, repetitions,pca_structure_correction, weighting, gene_conditioning, random_effects, pvalue_boot_threshold, regression_dataframe)
        }
    stopImplicitCluster()
    
    file_name = make_regression_file_name(snp_list, trim_type, condensing, random_effects, d_infer, repetitions, pca_structure_correction)
    if (write_table == "True"){
        # Write tables
        write.table(as.data.frame(results), file= file_name, quote=FALSE, sep='\t', col.names = NA)
    } else if (write_table != "True"){
        return(results) 
    }
}


