#!/bin/sh

# This script will make the cluster submission file for all insertion types
set -eu

SCRIPT=$1
PHENOTYPE=$2
PARTITION=$3
NCPU=$4

echo "bash $PWD/submission_scripts/$SCRIPT vd_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/vd_${PHENOTYPE}_submissions
bash $PWD/submission_scripts/$SCRIPT dj_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/dj_${PHENOTYPE}_submissions
bash $PWD/submission_scripts/$SCRIPT vj_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/vj_${PHENOTYPE}_submissions
bash $PWD/submission_scripts/$SCRIPT total_${PHENOTYPE} $PARTITION $NCPU > $PWD/submission_outputs/total_${PHENOTYPE}_submissions" > $PWD/submission_scripts/phenotype_functions/phenotype_class/submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}

echo "$PWD/submission_scripts/phenotype_functions/phenotype_class/submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}"
