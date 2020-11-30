#!/bin/sh

set -eu

PROJECT_PATH="/home/mrussel2"
OUTPUT_PATH="/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output"
TRIM_TYPE=$1
PCA=$2


get_path_command="Rscript $PROJECT_PATH/tcr-gwas/trimming_regression/scripts/get_regression_output_file_path.R $TRIM_TYPE $PCA $OUTPUT_PATH $PROJECT_PATH"
REGRESSION_OUTPUT_PATH=$($get_path_command)
echo $REGRESSION_OUTPUT_PATH
