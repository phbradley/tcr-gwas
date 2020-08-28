#!/bin/zsh

for SCRIPT in $(find $PWD -name "run_script_on_cluster.sh")
do
    cd $(dirname $SCRIPT)
    sbatch $SCRIPT
done
