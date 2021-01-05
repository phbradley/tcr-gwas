#!/bin/sh
set -eu

PHENOTYPE=$1

source $PWD/config/config.sh $PHENOTYPE
for JOB in {01..35481497..1000}; do
    mkdir -p $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}/${PHENOTYPE}_${JOB}
    cd $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}/${PHENOTYPE}_${JOB}
    echo "#!/bin/bash
    source $HOME/miniconda3/etc/profile.d/conda.sh
    conda activate tcr-gwas 
    set -eu
    Rscript $PROJECT_PATH/tcr-gwas/gwas_regressions/src/execute_regressions.R ${JOB} $PHENOTYPE $(echo '$1') 
    echo "done" > run.sentinel" > run_regressions_on_cluster.sh
    chmod +x run_regressions_on_cluster.sh
    chmod +x $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}/${PHENOTYPE}_${JOB}/run_regressions_on_cluster.sh
    cp $PROJECT_PATH/tcr-gwas/gwas_regressions/config/config.R $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}/${PHENOTYPE}_${JOB}/
    cd $PROJECT_PATH/tcr-gwas/gwas_regressions
done
echo "Submission directories created"
