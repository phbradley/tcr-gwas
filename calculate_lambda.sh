#!/bin/bash

# This script will continuously sumbit cluster jobs until all jobs are complete
set -eu

PHENOTYPE=$1
PHENOTYPE_CLASS=$2
PARTITION=$3
NCPU=$4

source $PWD/config/config.sh $PHENOTYPE

if [ "$PHENOTYPE_CLASS" = "trimming" ]; then
    BOOTSTRAP_OUTPUT_PATH="${OUTPUT_PATH}/results/bootstrap_lambda_analysis/${PHENOTYPE}_100"
    BOOTSTRAP_FILES=$(ls $BOOTSTRAP_OUTPUT_PATH | wc -l)
    if [ $BOOTSTRAP_FILES -lt 100 ]; then
        COMMAND="source $PROJECT_PATH/tcr-gwas/submission_scripts/submit_all_random_bootstraps_analyses.sh $PHENOTYPE $PARTITION $NCPU"
        $COMMAND
    fi
    while [ $BOOTSTRAP_FILES -lt 100 ]; do 
        sleep 5m 
        BOOTSTRAP_FILES=$(ls $BOOTSTRAP_OUTPUT_PATH | wc -l)
    done
    LAMBDA_CALC="Rscript $PROJECT_PATH/tcr-gwas/analysis_scripts/calculate_lambdas_gene_conditioned.R $PHENOTYPE $NCPU"
    echo $LAMBDA_CALC
    $LAMBDA_CALC
else
    LAMBDA_CALC="Rscript $PROJECT_PATH/tcr-gwas/analysis_scripts/calculate_lambdas.R $PHENOTYPE $NCPU"
    echo $LAMBDA_CALC
    $LAMBDA_CALC
fi
echo "done calculating lambdas for $PHENOTYPE"
