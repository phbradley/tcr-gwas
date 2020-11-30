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

bash $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/make_combination_directories.sh $TRIM_TYPE $PCA $PCA_TYPE

bash $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/continuous_submit.sh $TRIM_TYPE $PCA $PARTITION $CPU_COUNT

cd $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/${TRIM_TYPE}_${PCA}

sbatch -c $CPU_COUNT -p $PARTITION -q $PARTITION $PROJECT_PATH/tcr-gwas/trimming_regression/scripts/compile_regressions.sh $TRIM_TYPE $PCA $PCA_TYPE $PROJECT_PATH $OUTPUT_PATH > compile_job_id 

while [ $(squeue -u $USER -j $(cut -c21- compile_job_id) | wc -l) -eq 2 ]; do
    sleep 10m
    echo "still compiling regressions"
done 

echo "finished compiling regressions!"

rm -r $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/${TRIM_TYPE}_${PCA}
rm -r $REGRESSION_OUTPUT_PATH

