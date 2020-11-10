#!/bin/sh
set -eu

TRIM_TYPE=$1
PCA=$2
PCA_TYPE=$3

source parameters.sh $TRIM_TYPE $PCA

if [ ! -d "$PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories" ]; then
    mkdir $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories 
fi

if [ ! -d "$PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/${TRIM_TYPE}_${PCA}" ]; then
    mkdir $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/${TRIM_TYPE}_${PCA}
fi

for number in {01..35481497..1000}; do
    mkdir $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/${TRIM_TYPE}_${PCA}/${TRIM_TYPE}_${PCA}_${number}
    cd  $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/${TRIM_TYPE}_${PCA}/${TRIM_TYPE}_${PCA}_${number}
    echo "#!/bin/bash
    source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
    conda activate r
    set -eu
    Rscript $PROJECT_PATH/tcr-gwas/trimming_regression/scripts/run_cluster_random.R ${number} $TRIM_TYPE $PCA $(echo '$1') $PROJECT_PATH $OUTPUT_PATH $PCA_TYPE
    echo "done" > run.sentinel" > run_script_on_cluster.sh
    chmod +x run_script_on_cluster.sh
    chmod +x $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/${TRIM_TYPE}_${PCA}/${TRIM_TYPE}_${PCA}_${number}/run_script_on_cluster.sh
    cd $PROJECT_PATH/tcr-gwas/trimming_regression
done
echo "Submission directories created"
