library('doParallel')
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_regression_functions.R')

args = commandArgs(trailingOnly=TRUE)


trim_type = args[1]
random_effects = 'False'
condensing = 'by_patient'
d_infer = 'False'
repetitions = 0
pca_structure_correction = 'False'

data_files = find_files(trim_type, random_effects, condensing, d_infer, repetitions, pca_structure_correction)

registerDoParallel(cores=as.numeric(args[2]))
    results = foreach(file = data_files, .combine='rbind') %dopar% {
        RcppParallel::setThreadOptions(1L) 
        compile_all_data_from_cluster_for_parallel(file, data_files, trim_type, random_effects, condensing, d_infer, repetitions, pca_structure_correction)
    }
stopImplicitCluster()

subset_and_write_df(results, trim_type, random_effects, condensing, repetitions, pca_structure_correction, ncpus = as.numeric(args[2]))
