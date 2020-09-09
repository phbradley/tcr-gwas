#!/bin/zsh

set -eu

for CONFIG in $(find $PWD/cluster_job_directories/v_trim* -name "run_script_on_cluster.sh")
do
    cd $(dirname $CONFIG)
    if test -n "$(find . -maxdepth 1 -name '*sentinel' -print -quit)"
    then
        echo "\n\n*** Skipping $CONFIG because sentinel file found."
        continue
    fi
    if test -n "$(find . -maxdepth 1 -name 'slurm*.out' -print -quit)"
    then
        echo "\n\n*** Skipping $CONFIG because process is currently running."
        continue
    fi
    echo "\n*** Starting $CONFIG"
    COMMAND="sbatch -c 4 $CONFIG"
    echo "Going to run \`$COMMAND\`"
    $COMMAND
    echo "*** Completed $CONFIG"
done
