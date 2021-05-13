# Configuration Files
These files contain variables for repository path, output data path, data paths, etc. 
The variables must be changed to be computer/project specific. 

Specifically, the following two changes are required:

1. `PROJECT_PATH` (location of the `tcr-gwas` repository) and `OUTPUT_PATH` (location where cluster outputs should be stored) variables within [config.sh](config.sh) and [config.R](config.R) must be changed
2. All data paths should be changed within the [file_paths.R](file_paths.R), [validation_file_paths_alpha.R](validation_file_paths_alpha.R), and [validation_file_paths_beta.R](validation_file_paths_beta.R)

