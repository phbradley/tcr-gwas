#!/bin/bash

source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu

PHENOTYPE=$1
NCPU=$2

source $PWD/config/config.sh $PHENOTYPE

Rscript $PROJECT_PATH/tcr-gwas/gwas_regressions/scripts/compile_regressions.R $PHENOTYPE $NCPU
