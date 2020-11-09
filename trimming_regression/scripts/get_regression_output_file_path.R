args = commandArgs(trailingOnly=TRUE)
project_path = args[8]

source(paste0(project_path, '/tcr-gwas/trimming_regression/scripts/compile_data_functions.R'))
make_regression_file_path(trim_type = args[1], condensing = args[2], random_effects = args[3], d_infer = args[4], repetitions = args[5], pca_structure_correction = args[6], output_path = args[7])
