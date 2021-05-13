# Genome-wide analysis scripts

This directory contains files for running the analyses genome-wide for a variety of phenotypes: 

* [Analysis execution script](execute_regressions.R) which is sourced by top-level analysis run scripts (i.e. ../submit_phenotype.sh and ../submit_phenotype_class.sh) 
* [General analysis functions](regression_functions.R) which is sourced by most other scripts
* Cluster run compilation script in [R](compile_regressions.R) and [shell](compile_regressions.sh) 
* [Cluster job output file path generator](get_regression_output_path.R)

Phenotype specific functions and variables are located in the following [directory](phenotype_functions)


