#!/bin/sh
set -eu

SCRIPT=$1
PHENOTYPE=$2
PARTITION=$3
NCPU=$4

echo "bash $SCRIPT vd_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/vd_${PHENOTYPE}_submissions
bash $SCRIPT dj_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/dj_${PHENOTYPE}_submissions
bash $SCRIPT vj_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/vj_${PHENOTYPE}_submissions
bash $SCRIPT total_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/total_${PHENOTYPE}_submissions" > $PWD/submission_scripts/phenotype_functions/phenotype_class/submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}

echo "$PWD/submission_scripts/phenotype_functions/phenotype_class/submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}"
