#!/bin/zsh

set -eu

for CONFIG in $(find $PWD/cluster_job_directories/* -name "run_script_on_cluster.sh")
do
    cd $(dirname $CONFIG)
    cp /home/mrussel2/tcr-gwas/trimming_regression/run_script_on_cluster.sh .
    chmod +x $CONFIG
done
