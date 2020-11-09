#!/bin/sh
set -eu

TRIM_TYPE=$1
PCA=$2

source parameters.sh $TRIM_TYPE $PCA

if [ ! -d "$PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories" ]; then
    mkdir $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories 
fi

if [ ! -d "$PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/$TRIM_TYPE" ]; then
    mkdir $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/$TRIM_TYPE
fi

for number in {01..35481497..1000}; do
    mkdir $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/$TRIM_TYPE/${TRIM_TYPE}_${number}
    cd  $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/$TRIM_TYPE/${TRIM_TYPE}_${number}
    echo "#!/bin/bash
    source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
    conda activate r
    set -eu
    Rscript $PROJECT_PATH/tcr-gwas/trimming_regression/scripts/run_cluster_random.R ${number} $TRIM_TYPE $(echo '$2') $(echo '$1') $PROJECT_PATH $OUTPUT_PATH
    echo "done" > run.sentinel" > run_script_on_cluster.sh
    chmod +x run_script_on_cluster.sh
    chmod +x $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/$TRIM_TYPE/${TRIM_TYPE}_${number}/run_script_on_cluster.sh
    cd $PROJECT_PATH/tcr-gwas/trimming_regression
done
echo "Submission directories created"
