#!/bin/sh

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate tcr-gwas 
set -eu

PHENOTYPE=$1
NCPU=$2
CHAIN=$3

Rscript $PWD/analysis_scripts/validation_cohort_analysis/execute_validation_data_regressions.R $PHENOTYPE $NCPU $CHAIN
