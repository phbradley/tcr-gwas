library('foreach')
library('doParallel')
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_regression_functions.R')

args = commandArgs(trailingOnly=TRUE)

stopifnot(args[1] %in% c('trimming', 'insertion')
types = ifelse(args[1] == 'trimming', c('v_trim', 'd0_trim', 'd1_trim', 'j_trim'), c('vd_insert', 'vj_insert', 'dj_insert'))

# Use this, for running different trimming type compiles in parallel, but each
# is indidividually run sequentially

registerDoParallel(cores=args[2])
    foreach(trim_type=types) %dopar% {
        RcppParallel::setThreadOptions(1L) 
        compile_all_data_from_cluster_sequential(trim_type = trim_type, pca_structure_correction = 'True') }
stopImplicitCluster()
