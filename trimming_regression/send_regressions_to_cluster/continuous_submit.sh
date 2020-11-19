#!/bin/sh

set -eu

TRIM_TYPE=$1
PARTITION=$3
CPU_COUNT=$4
PCA=$2
source parameters.sh $TRIM_TYPE $PCA    

TOTAL_COUNT=$(ls $REGRESSION_OUTPUT_PATH | wc -l)
#EXPECTED_TOTAL=$(ls $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/$TRIM_TYPE | wc -l)
EXPECTED_TOTAL=35479
while [ $TOTAL_COUNT -lt $EXPECTED_TOTAL ]; do
    JOB_COUNT=$(squeue -u $USER -p $PARTITION | wc -l)
    echo "There are currently $JOB_COUNT jobs submitted to the normal cluster for $TRIM_TYPE"
    while [ $JOB_COUNT -le 1000 ]; do
        TOTAL_COUNT=$(ls $REGRESSION_OUTPUT_PATH | wc -l)
        if [[ "$TOTAL_COUNT" == "$EXPECTED_TOTAL" ]]; then
            break
        fi
        for CONFIG in $(find $PROJECT_PATH/tcr-gwas/trimming_regression/send_regressions_to_cluster/cluster_directories/${TRIM_TYPE}_${PCA}/ -name "run_script_on_cluster.sh" -exec ls {} + ); do
            cd $(dirname $CONFIG)
            if test -n "$(find . -maxdepth 1 -name '*sentinel' -print -quit)"
            then
                echo "\n\n*** Skipping $CONFIG because sentinel file found."
                continue
            fi
            if test -n "$(find . -maxdepth 1 -name 'job_id_running' -print -quit)"
            then
    	        if (($(squeue -u $USER -j $(cut -c21- job_id_running) | wc -l) == 2 ))
    	        then 
    		    echo "\n\n*** Skipping $CONFIG because process is currently running."
            	    continue
    	        fi
            fi
            COMMAND="sbatch -c $CPU_COUNT -p $PARTITION -q $PARTITION $CONFIG $CPU_COUNT"
            $COMMAND > job_id_running
            echo "Running \`$COMMAND\`"
            JOB_COUNT=$(squeue -u $USER -p $PARTITION | wc -l)
            if (( $JOB_COUNT >= 1000 ))
            then 
                break
            fi
        done
    done
    echo "There are currently $JOB_COUNT jobs submitted to the normal cluster and $TOTAL_COUNT jobs completed overall for $TRIM_TYPE"
    TOTAL_COUNT=$(ls $REGRESSION_OUTPUT_PATH | wc -l)
done
echo "Finished with regressions"
