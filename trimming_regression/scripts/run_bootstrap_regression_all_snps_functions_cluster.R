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
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_regression_data_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/scripts/execute_regression_function.R")

# This function runs regressions on the cluster
run_snps_trimming_snp_list_cluster <- function(snp_list, genotype_list, trim_type, gene_type, condensing, gene_conditioning, weighting, random_effects, repetitions, pca_structure_correction, write_table, ncpus, d_infer, maf_cutoff, data_file_path){
    regression_dataframe = data.frame()

    # import condensed trimming file
    if (condensing == 'by_patient'){
        if (d_infer == 'False'){
            assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/by_patient_condensed_data_all_patients_NO_d_infer.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        } else {
            assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/by_patient_condensed_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        }
        names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"
    } else if (condensing == 'gene_cross'){
        trimming_data = compile_trimming_data_cross()
        trimming_data = trimming_data %>% filter(gene_class == gene_type)
        names(trimming_data)[names(trimming_data) == 'weighted_gene_count'] <- paste0('weighted_', gene_type, '_count')
        names(trimming_data)[names(trimming_data) == 'gene'] <- paste0(gene_type)
    } else {
        if (d_infer == 'False'){
            assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/condensed_", trim_type, "_data_all_patients_NO_d_gene_infer.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        } else {
            assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/condensed_", trim_type, "_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        }
        names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"
    }

    #filter out snps with maf below cutoff
    if (maf_cutoff == "False"){
        list_of_snps = as.numeric(colnames(genotype_list)[colnames(genotype_list) != 'localID'])
    } else {
        maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]

        maf_data_filtered = maf_data %>% filter(maf >= maf_cutoff)

        maf_data_snps = gsub('snp', '', maf_data_filtered$snp)

        list_of_snps = as.numeric(intersect(colnames(genotype_list), maf_data_snps))
    }
    

    #results = do.call(rbind, mclapply(as.numeric(colnames(genotype_list)), execute_regression, snp_list, genotype_list, trim_type, gene_type, condensing, trimming_data, repetitions, weighting, gene_conditioning, random_effects, regression_dataframe, mc.cores = ncpus))
    count = 0
    registerDoParallel(cores=ncpus)
    results = foreach(snp=list_of_snps, .combine='rbind') %dopar% {
        RcppParallel::setThreadOptions(1L) 
        execute_regression(snp, list_of_snps, snp_list, genotype_list, trim_type, gene_type, condensing, trimming_data, repetitions,pca_structure_correction, weighting, gene_conditioning, random_effects, regression_dataframe)
        }
    stopImplicitCluster()
    
    if (random_effects == 'True'){
        if (d_infer == 'False') {
            file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/NO_d_infer/',trim_type, '/', repetitions, '_bootstraps', trim_type, '_',snp_list$snp[1], '-', snp_list$snp[nrow(snp_list)],'_snps_regression_with_weighting_condensing_', condensing, '_with_random_effects_', repetitions, '_bootstraps_NO_d_infer')
        } else {
            file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/',trim_type, '/', repetitions, '_bootstraps/', trim_type, '_',snp_list$snp[1], '-', snp_list$snp[nrow(snp_list)],'_snps_regression_with_weighting_condensing_', condensing, '_with_random_effects_', repetitions, '_bootstraps')
        }
    } else {
        if (d_infer == 'False') {
            file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/NO_d_infer/',trim_type, '/', repetitions, '_bootstraps/', trim_type, '_',snp_list$snp[1], '-', snp_list$snp[nrow(snp_list)],'_snps_regression_with_weighting_condensing_', condensing, '_NO_random_effects_', repetitions, '_bootstraps_NO_d_infer')
        } else {
            file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/',trim_type, '/', repetitions, '_bootstraps/', trim_type, '_',snp_list$snp[1], '-', snp_list$snp[nrow(snp_list)],'_snps_regression_with_weighting_condensing_', condensing, '_NO_random_effects_', repetitions, '_bootstraps')
        }
    }

    if (pca_structure_correction == 'True'){
        file_name = paste0(file_name, '_WITH_pca_structure_correction.tsv')
    }else{
        file_name = paste0(file_name, '.tsv')
    }
            
    if (write_table == "True"){
        # Write tables
        write.table(as.data.frame(results), file= file_name, quote=FALSE, sep='\t', col.names = NA)
    } else if (write_table != "True"){
        return(results) 
    }
}


