#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL
source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu

TRIM_TYPE=$1
PCA_STRUCTURE_CORRECTION=$2
PCA_TYPE=$3
PROJECT_PATH=$4
OUTPUT_PATH=$5

Rscript $PROJECT_PATH/tcr-gwas/trimming_regression/scripts/compile_regressions_together.R $TRIM_TYPE $PCA_STRUCTURE_CORRECTION $PCA_TYPE $PROJECT_PATH $OUTPUT_PATH
