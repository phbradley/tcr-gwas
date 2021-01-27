#!/bin/sh
set -eu

SCRIPT=$1
PHENOTYPE=$2
PARTITION=$3
NCPU=$4

echo "bash $PWD/submission_scripts/$SCRIPT v_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/v_${PHENOTYPE}_submissions
bash $PWD/submission_scripts/$SCRIPT j_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/j_${PHENOTYPE}_submissions
bash $PWD/submission_scripts/$SCRIPT d0_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/d0_${PHENOTYPE}_submissions
bash $PWD/submission_scripts/$SCRIPT d1_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/d1_${PHENOTYPE}_submissions" > $PWD/submission_scripts/phenotype_functions/phenotype_class/submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}

echo "${PWD}/submission_scripts/phenotype_functions/phenotype_class/submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}"
