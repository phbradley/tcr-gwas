#!/bin/sh

set -eu

PROJECT_PATH="/home/mrussel2"
OUTPUT_PATH="/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output"
TRIM_TYPE=$1
PCA=$2

TRIMMING="v_trim j_trim d0_trim d1_trim"
INSERTION="vd_insert dj_insert vj_insert"

if [[ $TRIMMING == *$TRIM_TYPE* ]]
then
    CONDENSE="by_gene"
    RANDOM_EFFECT="True"
    D_INFER="True"
    REPS="100"
elif [[ $INSERTION == *$TRIM_TYPE* ]]
then
    CONDENSE="by_patient"
    RANDOM_EFFECT="False"
    D_INFER="False"
    REPS="0"
else 
    echo "Trim type entered is not valide"
fi

get_path_command="Rscript $PROJECT_PATH/tcr-gwas/trimming_regression/scripts/get_regression_output_file_path.R $TRIM_TYPE $CONDENSE $RANDOM_EFFECT $D_INFER $REPS $PCA $OUTPUT_PATH $PROJECT_PATH"
REGRESSION_OUTPUT_PATH=$($get_path_command)
echo $REGRESSION_OUTPUT_PATH
