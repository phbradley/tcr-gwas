#!/bin/sh
set -eu

SCRIPT=$1
PHENOTYPE="tcr_div"
PARTITION=$2
NCPU=$3

bash $PWD/submission_scripts/make_submission_file_single.sh $SCRIPT $PHENOTYPE $PARTITION $NCPU
