#!/bin/sh
set -eu

SCRIPT=$1
PHENOTYPE=$2
PARTITION=$3
NCPU=$4

echo "bash $SCRIPT v_${PHENOTYPE} $PARTITION $NCPU > submission_outputs/v_${PHENOTYPE}_submissions
bash $SCRIPT j_${PHENOTYPE} $PARTITION $NCPU > submission_outputs/j_${PHENOTYPE}_submissions
bash $SCRIPT d0_${PHENOTYPE} $PARTITION $NCPU > submission_outputs/d0_${PHENOTYPE}_submissions
bash $SCRIPT d1_${PHENOTYPE} $PARTITION $NCPU > submission_outputs/d1_${PHENOTYPE}_submissions" > submit_${SCRIPT}_${PHENOTYPE}_${PARTITION}


