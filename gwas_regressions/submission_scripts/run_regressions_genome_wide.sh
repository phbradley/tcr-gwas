#!/bin/sh

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate tcr-gwas 
set -eu

PHENOTYPE=$1
PARTITION=$2
NCPU=$3

source $PWD/config/config.sh $PHENOTYPE

bash $PROJECT_PATH/tcr-gwas/gwas_regressions/submission_scripts/make_job_directories.sh $PHENOTYPE

bash $PROJECT_PATH/tcr-gwas/gwas_regressions/submission_scripts/submit_cluster_jobs_continuously.sh $PHENOTYPE $PARTITION $NCPU

cp $PWD/config/config.sh $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}
cp $PWD/config/config.R $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}
cd $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}

sbatch -c $NCPU -p $PARTITION -q $PARTITION $PROJECT_PATH/tcr-gwas/gwas_regressions/scripts/compile_regressions.sh $PHENOTYPE $NCPU > compile_job_id

while [ $(squeue -u $USER -j $(cut -c21- compile_job_id) | wc -l) -eq 2 ]; do
    echo "still compiling regression files"
    sleep 20m
done

echo "finished compiling regression files!"

rm -r $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}
rm -r $REGRESSION_OUTPUT_PATH