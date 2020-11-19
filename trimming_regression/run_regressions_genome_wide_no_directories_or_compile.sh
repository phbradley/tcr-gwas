#!/bin/sh

source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu

TRIM_TYPE=$1
PCA=$2
PARTITION=$3
CPU_COUNT=$4
PCA_TYPE=$5

source parameters.sh $TRIM_TYPE $PCA

bash $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/continuous_submit.sh $TRIM_TYPE $PCA $PARTITION $CPU_COUNT
echo "finished regressions, need to run the compile script"
