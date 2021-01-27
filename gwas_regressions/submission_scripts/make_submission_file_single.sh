#!/bin/sh
set -eu

SCRIPT=$1
PHENOTYPE=$2
PARTITION=$3
NCPU=$4

echo "bash $PWD/submission_scripts/$SCRIPT ${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/${PHENOTYPE}_submissions" > $PWD/submission_scripts/phenotype_functions/phenotype_class/submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}

echo "${PWD}/submission_scripts/phenotype_functions/phenotype_class/submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}"
