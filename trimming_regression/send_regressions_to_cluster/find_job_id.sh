#!/bin/zsh

set -eu

for CONFIG in $(find $PWD/random_effects_cluster_job_directory_100_boots_subset/$1/$1* -name "job_id_running"); do
    cd $(dirname $CONFIG)
    if (($(cut -c21- job_id_running) == $2))
	then 
	    echo "$(dirname $CONFIG)"
        break
	fi
done
