library('foreach')
library('doParallel')
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_regression_functions.R')

# Use this, for running different trimming type compiles in parallel, but each
# is indidividually run sequentially

registerDoParallel(cores=3)
    foreach(trim_type=c('vd_insert', 'vj_insert', 'dj_insert')) %dopar% { 
    #foreach(trim_type=c('v_trim', 'd0_trim', 'd1_trim', 'j_trim')) %dopar% {
        RcppParallel::setThreadOptions(1L) 
        compile_all_data_from_cluster_sequential(trim_type = trim_type, random_effects = 'False', condensing = 'by_patient', d_infer = 'False', repetitions = 100, pca_structure_correction = 'True') }
stopImplicitCluster()
