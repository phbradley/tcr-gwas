#!/bin/bash
set -eu

SCRIPT=$1
PHENOTYPE="trim_naive"
PARTITION=$2
NCPU=$3

bash $PWD/submission_scripts/make_submission_file.sh $SCRIPT $PHENOTYPE $PARTITION $NCPU
