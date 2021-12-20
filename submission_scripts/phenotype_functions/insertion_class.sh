#!/bin/bash
set -eu

SCRIPT=$1
PHENOTYPE="insert"
PARTITION=$2
NCPU=$3

bash $PWD/submission_scripts/make_submission_file_inserts.sh $SCRIPT $PHENOTYPE $PARTITION $NCPU
