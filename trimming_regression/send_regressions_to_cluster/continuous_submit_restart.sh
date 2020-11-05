#!/bin/zsh

set -eu

TOTAL_COUNT=$(ls /fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/NO_d_infer/$1/0_bootstraps | wc -l)
while (( $TOTAL_COUNT < 35500 ))
do
    JOB_COUNT=$(squeue -u $USER -p restart-new | wc -l)
    echo "There are currently $JOB_COUNT jobs submitted to the restart cluster for $1"
    while (( $JOB_COUNT < 1000 ))
    do
        for CONFIG in $(find /home/mrussel2/tcr-gwas/trimming_regression/send_regressions_to_cluster/random_effects_cluster_job_directory/$1/ -name "run_script_on_cluster.sh" -exec ls {} + ); do
            cd $(dirname $CONFIG)
            if test -n "$(find . -maxdepth 1 -name '*sentinel' -print -quit)"
            then
                echo "\n\n*** Skipping $CONFIG because sentinel file found."
                continue
            fi
            if test -n "$(find . -maxdepth 1 -name 'job_id_running' -print -quit)"
            then
    	        if (($(squeue -u $USER -j $(cut -c21- job_id_running) | wc -l) == 2))
    	        then 
    		    echo "\n\n*** Skipping $CONFIG because process is currently running."
            	    continue
    	        fi
            fi
            COMMAND="sbatch -p restart-new -c $2 $CONFIG $2"
            $COMMAND > job_id_running
            echo "Running \`$COMMAND\`"
            JOB_COUNT=$(squeue -u $USER -p restart-new | wc -l)
            if (( $JOB_COUNT >= 1000 ))
            then 
                break
            fi
        done
    done
    echo "There are currently $JOB_COUNT jobs submitted to the restart cluster and $TOTAL_COUNT jobs completed overall for $1"
    sleep 10m
    TOTAL_COUNT=$(ls /fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/NO_d_infer/$1/0_bootstraps | wc -l)
done
echo "ALL DONE!!!"
