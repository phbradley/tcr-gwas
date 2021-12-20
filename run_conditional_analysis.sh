#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate tcr-gwas 
set -eu

GENE=$1
PHENOTYPE=$2
NCPU=$3

Rscript $PWD/analysis_scripts/conditional_analysis_scripts/execute_conditional_analysis.R $GENE $PHENOTYPE $NCPU 
