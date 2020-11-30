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


source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts_by_race/regression_functions.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts_by_race/bootstrap_functions.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts_by_race/compile_data_functions.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts_by_race/execute_regression_function.R"))

# This function runs regressions on the cluster
run_snps_trimming_snp_list_cluster_race_covariate <- function(snp_list, genotype_list, trim_type, pca_structure_correction, pca_type, by_race, write_table){
    
    stopifnot(trim_type %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'vj_insert'))
    regression_dataframe = data.frame()

    # import condensed trimming file
    condensed_trimming_data = compile_condensed_trimming_data(trim_type, CONDENSING)
    
    # filter out small repertoires
    trimming_data = remove_small_repertoire_observations(condensed_trimming_data, 
                                                         productive_log10_count_cutoff= 4.25, 
                                                         NOT_productive_log10_count_cutoff=3.5, 
                                                         trim_type) 

    if (by_race == 'True'){
        ethnicity_data = fread(file = '/home/mrussel2/tcr-gwas/_ignore/race_pcs_18Nov2020.txt')   
        trimming_data = merge(trimming_data, ethnicity_data[,-c('scanID')], by = 'localID')  
    }

    #filter out snps with maf below cutoff
    list_of_snps = filtered_snps_by_maf(MAF_CUTOFF, genotype_list)
    
    print('All data compiled. Starting regressions.')

    count = 0
    registerDoParallel(cores=NCPU)
    results = foreach(snp=list_of_snps, .combine='rbind') %dopar% {
        RcppParallel::setThreadOptions(1L) 
        execute_regression(snp, 
                           list_of_snps, 
                           snp_list, 
                           genotype_list, 
                           trim_type, 
                           trimming_data, 
                           pca_structure_correction, 
                           pca_type,
                           by_race,
                           regression_dataframe)
        }
    stopImplicitCluster()
    
    file_name = make_regression_file_name(snp_list, 
                                          trim_type, 
                                          pca_structure_correction, 
                                          by_race)

    if (write_table == "True"){
        write.table(as.data.frame(results), file= file_name, quote=FALSE, sep='\t', col.names = NA)
    } else if (write_table != "True"){
        return(results) 
    }
}


