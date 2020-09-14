#!/bin/zsh

set -eu

for number in {1..9..1}; do for CONFIG in $(find $PWD/random_effects_cluster_job_directory/$1/$1_$number* -name "run_script_on_cluster.sh"); do
    cd $(dirname $CONFIG)
    if test -n "$(find . -maxdepth 1 -name '*sentinel' -print -quit)"
    then
        echo "\n\n*** Skipping $CONFIG because sentinel file found."
        continue
    fi
    if test -n "$(find . -maxdepth 1 -name 'job_id_running' -print -quit)"
    then
	    if (($(sacct -j $(cut -c21- job_id_running) | wc -l) > 2 && $(sacct -j $(cut -c21- job_id_running) | wc -l) < 5))
	    then 
		echo "\n\n*** Skipping $CONFIG because process is currently running."
        	continue
	    fi
    fi
    COMMAND="sbatch -c $2 -p restart-new $CONFIG $2"
    $COMMAND > job_id_running
    echo "Running \`$COMMAND\`"
done
done
