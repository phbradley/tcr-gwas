#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu

PHENOTYPE=$1
NCPU=$2

source $PWD/config/config.sh $PHENOTYPE

Rscript $PROJECT_PATH/tcr-gwas/scripts/compile_regressions.R $PHENOTYPE $NCPU
