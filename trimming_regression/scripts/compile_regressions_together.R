args = commandArgs(trailingOnly=TRUE)

project_path = args[3]
output_path = args[4]

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_regression_functions.R')

# Use this, for running different trimming type compiles in parallel, but each
# is indidividually run sequentially

compile_all_data_from_cluster_sequential(trim_type = args[1], pca_structure_correction = args[2]) 
