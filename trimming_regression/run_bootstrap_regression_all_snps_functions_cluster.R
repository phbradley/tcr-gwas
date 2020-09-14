library('plyr')
library("tidyverse")
library('foreach')
library('doParallel')
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)


source("/home/mrussel2/tcr-gwas/trimming_regression/simple_trimming_regression_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/lmer_trimming_regression_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/bootstrap_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/compile_regression_data_functions.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/execute_regression_function.R")


run_snps_trimming_snp_list_cluster <- function(snp_list, genotype_list, trim_type, gene_type, condensing, gene_conditioning, weighting, random_effects, repetitions, write_table, ncpus){
    regression_dataframe = data.frame()

    # import condensed trimming file
    if (condensing == 'by_patient'){
        assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/by_patient_condensed_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"
    } else if (condensing == 'gene_cross'){
        trimming_data = compile_trimming_data_cross()
        trimming_data = trimming_data %>% filter(gene_class == gene_type)
        names(trimming_data)[names(trimming_data) == 'weighted_gene_count'] <- paste0('weighted_', gene_type, '_count')
        names(trimming_data)[names(trimming_data) == 'gene'] <- paste0(gene_type)
    } else if (condensing == 'phil'){
        assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/by_patient_condensed_data_all_patients_phil.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
    } else {
        assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/condensed_", trim_type, "_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"
    }

    #results = do.call(rbind, mclapply(as.numeric(colnames(genotype_list)), execute_regression, snp_list, genotype_list, trim_type, gene_type, condensing, trimming_data, repetitions, weighting, gene_conditioning, random_effects, regression_dataframe, mc.cores = ncpus))
    registerDoParallel(cores=ncpus)
    results = foreach(snp=as.numeric(colnames(genotype_list)), .combine='rbind') %dopar% {
        RcppParallel::setThreadOptions(1L) 
        execute_regression(snp = snp, snp_list, genotype_list, trim_type, gene_type, condensing, trimming_data, repetitions, weighting, gene_conditioning, random_effects, regression_dataframe)
        }
    stopImplicitCluster()
    
    if (random_effects == 'True'){
         file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/',trim_type, '/', trim_type, '_',snp_list$snp[1], '-', snp_list$snp[nrow(snp_list)],'_snps_regression_with_weighting_condensing_by_gene_with_random_effects_', repetitions, '_bootstraps.tsv')
    } else {
         file_name = paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/',trim_type, '/', trim_type, '_',snp_list$snp[1], '-', snp_list$snp[nrow(snp_list)],'_snps_regression_with_weighting_condensing_by_gene_NO_random_effects_', repetitions, '_bootstraps.tsv')
    }
            
    if (write_table == "True"){
        # Write tables
        write.table(as.data.frame(results), file= file_name, quote=FALSE, sep='\t', col.names = NA)
    } else if (write_table != "True"){
        return(results)
    }
}


