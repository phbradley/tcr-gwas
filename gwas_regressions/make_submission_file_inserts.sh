#!/bin/sh
set -eu

SCRIPT=$1
PHENOTYPE=$2
PARTITION=$3
NCPU=$4

echo "bash $SCRIPT vd_${PHENOTYPE} $PARTITION $NCPU > submission_outputs/vd_${PHENOTYPE}_submissions
bash $SCRIPT dj_${PHENOTYPE} $PARTITION $NCPU > submission_outputs/dj_${PHENOTYPE}_submissions
bash $SCRIPT vj_${PHENOTYPE} $PARTITION $NCPU > submission_outputs/vj_${PHENOTYPE}_submissions
bash $SCRIPT total_${PHENOTYPE} $PARTITION $NCPU > submission_outputs/total_${PHENOTYPE}_submissions" > submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}


