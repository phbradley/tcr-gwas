#!/bin/bash

set -eu

PHENOTYPE=$1
PARTITION=$2
NCPU=$3

source $PWD/config/config.sh $PHENOTYPE

SCRIPT=$PROJECT_PATH/tcr-gwas/submission_scripts/bootstrap_analysis.sh 

COMMAND="sbatch -c $NCPU -p $PARTITION -q $PARTITION $SCRIPT $PHENOTYPE $NCPU"

COUNT=0
EXPECTED_COUNT=100

while [ $COUNT -lt $EXPECTED_COUNT ]; do
    $COMMAND
    COUNT=$[$COUNT+1]
done
