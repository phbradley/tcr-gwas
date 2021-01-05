#!/bin/sh

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate tcr-gwas
set -eu

PHENOTYPE_CLASS=$1
PARTITION=$2
NCPU=$3
default_script="run_regressions_genome_wide.sh"
SCRIPT=${4:-$default_script}

command=$(bash $PWD/submission_scripts/phenotype_functions/${PHENOTYPE_CLASS}_class.sh $SCRIPT $PARTITION $NCPU)
lines=$(cat $command | wc -l)

parallel --jobs $lines < $command

rm $command

