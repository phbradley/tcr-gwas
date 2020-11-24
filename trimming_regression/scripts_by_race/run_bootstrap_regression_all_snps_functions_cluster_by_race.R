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


source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/regression_functions.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/bootstrap_functions.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/execute_regression_function.R"))

# This function runs regressions on the cluster
run_snps_trimming_snp_list_cluster_by_race <- function(snp_list, genotype_list, trim_type, pca_structure_correction, pca_type, write_table){
    
    stopifnot(trim_type %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'vj_insert'))
    regression_dataframe = data.frame()

    # import condensed trimming file
    condensed_trimming_data = compile_condensed_trimming_data(trim_type, CONDENSING)
    
    # filter out small repertoires
    trimming_data = remove_small_repertoire_observations(condensed_trimming_data, 
                                                         productive_log10_count_cutoff= 4.25, 
                                                         NOT_productive_log10_count_cutoff=3.5, 
                                                         trim_type) 
    
    ethnicity_data = fread(file = '/home/mrussel2/tcr-gwas/_ignore/snp_data/ethnicity_data.csv')   

    together = merge(trimming_data, ethnicity_data, by.x = 'localID', by.y = 'id')


    #filter out snps with maf below cutoff
    list_of_snps = filtered_snps_by_maf(MAF_CUTOFF, genotype_list)
    
    print('All data compiled. Starting regressions.')
    
    together_regressions = data.frame()
    count = 0
    for (race_index in unique(together$race)){
        results = NULL
        if (nrow(together[race == race_index])>9){
            print(paste0('processing for ', race_index))
            registerDoParallel(cores=NCPU)
            results = foreach(snp=list_of_snps, .combine='rbind') %dopar% {
                RcppParallel::setThreadOptions(1L) 
                execute_regression(snp, 
                                   list_of_snps, 
                                   snp_list, 
                                   genotype_list[row.names(genotype_list) %in% together[race == race_index]$localID,], 
                                   trim_type, 
                                   together[race == race_index], 
                                   pca_structure_correction, 
                                   pca_type, 
                                   regression_dataframe)
            }
            stopImplicitCluster()
            results$race = race_index
            together_regressions = rbind(together_regressions, results)
        } else {
            print(paste0('skipping ', race_index))
        }
    }
    return(together_regressions)
}
