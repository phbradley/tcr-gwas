#!/bin/sh

source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu

PHENOTYPE=$1
PARTITION=$2
NCPU=$3

source config.sh $PHENOTYPE

cp config.sh $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}
cp config.R $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}
cd $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}

sbatch -c $NCPU -p $PARTITION -q $PARTITION $PROJECT_PATH/tcr-gwas/gwas_regressions/src/compile_regressions.sh $PHENOTYPE $NCPU > compile_job_id

while [ $(squeue -u $USER -j $(cut -c21- compile_job_id) | wc -l) -eq 2 ]; do
    echo "still compiling regression files"
    sleep 20m
done

echo "finished compiling regression files!"

rm -r $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}
rm -r $REGRESSION_OUTPUT_PATH
