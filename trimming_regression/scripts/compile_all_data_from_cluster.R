library('foreach')
library('doParallel')
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/manha_visualization.R')

#registerDoParallel(cores=4)
#    foreach(trim_type=c('d0_trim', 'd1_trim', 'j_trim', 'v_trim')) %dopar% {
#        RcppParallel::setThreadOptions(1L) 
#        compile_all_data_from_cluster(trim_type = trim_type, random_effects = 'True', bootstrap_count = 0) }
#stopImplicitCluster()

registerDoParallel(cores=1)
    foreach(trim_type=c('j_trim')) %dopar% {
        RcppParallel::setThreadOptions(1L) 
        compile_all_data_from_cluster(trim_type = trim_type, random_effects = 'True', bootstrap_count = 100) }
stopImplicitCluster()