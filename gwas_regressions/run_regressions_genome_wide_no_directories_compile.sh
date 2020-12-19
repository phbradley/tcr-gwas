#!/bin/sh

source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu

PHENOTYPE=$1
PARTITION=$2
NCPU=$3

source config.sh $PHENOTYPE

bash $PROJECT_PATH/tcr-gwas/gwas_regressions/src/submit_cluster_jobs_continuously.sh $PHENOTYPE $PARTITION $NCPU

