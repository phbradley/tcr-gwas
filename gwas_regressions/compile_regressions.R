source(paste0(getwd(), '/config.R'))

library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(vroom)
#TODO uncomment these
# args = commandArgs(trailingOnly=TRUE)

# PHENOTYPE <<- args[1]
# NCPU <<- args[2]

PHENOTYPE = 'd0_trim'
NCPU = 1

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/regression_functions.R'))

compile_regressions()
