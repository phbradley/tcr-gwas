#!/bin/sh

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate tcr-gwas
set -eu

PHENOTYPE=$1
PARTITION=$2
NCPU=$3
default_script="run_regressions_genome_wide.sh"
SCRIPT=${4:-$default_script}

bash $PWD/submission_scripts/$SCRIPT $PHENOTYPE $PARTITION $NCPU 

