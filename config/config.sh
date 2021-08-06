#!/bin/sh

set -eu

PROJECT_PATH=$HOME
OUTPUT_PATH="/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output"
PHENOTYPE=$1

get_regression_path="Rscript $PROJECT_PATH/tcr-gwas/scripts/get_regression_output_path.R $PHENOTYPE"

REGRESSION_OUTPUT_PATH=$($get_regression_path)

echo $REGRESSION_OUTPUT_PATH
