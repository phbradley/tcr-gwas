# Configuration Files
These files contain variables for repository path, output data path, data paths, etc. 
The variables must be changed to be computer/project specific. 

Specifically, the following two changes are required:

1. `PROJECT_PATH` (location of the `tcr-gwas` repository) and `OUTPUT_PATH` (location where cluster outputs should be stored) variables within `gwas-regressions/config/config.sh` and `gwas-regressions/config/config.R` must be changed
2. All data paths should be changed within the `gwas-regressions/config/file_paths.R` file

