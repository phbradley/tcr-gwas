#!/bin/sh

source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu

TRIM_TYPE=$1
PCA=$2
PARTITION=$3
CPU_COUNT=$4

source parameters.sh $TRIM_TYPE $PCA

bash $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/make_combination_directories.sh $TRIM_TYPE $PCA

bash $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/continuous_submit.sh $TRIM_TYPE $PCA $PARTITION $CPU_COUNT

Rscript $PROJECT_PATH/tcr-gwas/trimming_regression/scripts/compile_regressions_together.R $TRIM_TYPE $PCA $PROJECT_PATH $OUTPUT_PATH

rm -r $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories

