#!/bin/zsh

set -eu

for CONFIG in $(find $PWD/$3/$1* -name "job_id_running"); do
    cd $(dirname $CONFIG)
    if (($(cut -c21- job_id_running) == $2))
	then 
	    echo "$(dirname $CONFIG)"
        break
	fi
done
