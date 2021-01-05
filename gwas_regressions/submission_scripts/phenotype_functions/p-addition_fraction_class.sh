#!/bin/sh
set -eu

SCRIPT=$1
PHENOTYPE="pnucs_fraction"
PARTITION=$2
NCPU=$3

bash $PWD/submission_scripts/make_submission_file.sh $SCRIPT $PHENOTYPE $PARTITION $NCPU
