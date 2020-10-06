#!/bin/zsh

set -eu

for CONFIG in $(find $PWD/random_effects_cluster_job_directory/$1/$1_$2* -name "job_id_running"); do
    cd $(dirname $CONFIG)
    if (($(cut -c21- job_id_running) == $3))
	then 
	    echo "$(dirname $CONFIG)"
        break
	fi
done
