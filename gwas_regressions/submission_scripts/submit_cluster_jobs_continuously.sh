#!/bin/sh

set -eu

PHENOTYPE=$1
PARTITION=$2
NCPU=$3

source $PWD/config/config.sh $PHENOTYPE

FINISHED_REGRESSION_COUNT=$(ls $REGRESSION_OUTPUT_PATH | wc -l)
EXPECTED_REGRESSION_COUNT=35479
while [ $FINISHED_REGRESSION_COUNT -lt $EXPECTED_REGRESSION_COUNT ]; do
    RUNNING_JOBS_COUNT=$(squeue -u $USER -p $PARTITION | wc -l)
    # echo "There are currently $RUNNING_JOBS_COUNT jobs submitted to the $PARTITION cluster partition for $PHENOTYPE"
    while [ $RUNNING_JOBS_COUNT -le 1000 ]; do
        FINISHED_REGRESSION_COUNT=$(ls $REGRESSION_OUTPUT_PATH | wc -l)
        if [[ "$FINISHED_REGRESSION_COUNT" -ge "$EXPECTED_REGRESSION_COUNT" ]]; then
            break
        fi
        for JOB in $(find $PROJECT_PATH/tcr-gwas/gwas_regressions/current_cluster_job_directories/${PHENOTYPE}/ -name "run_regressions_on_cluster.sh" -exec ls {} + ); do
            cd $(dirname $JOB)
            if test -n "$(find . -maxdepth 1 -name '*sentinel' -print -quit)"
            then 
                # echo "\n\n*** Skipping $JOB because sentinel file found."
                continue
            fi
            if test -n "$(find . -maxdepth 1 -name 'job_id_running' -print -quit)"
            then 
                if (( $(squeue -u $USER -j $(cut -c21- job_id_running) | wc -l) == 2 ))
                then 
                    # echo "\n\n*** Skipping $JOB because process is currently running."
                    continue
                fi
            fi
            COMMAND="sbatch -c $NCPU -p $PARTITION -q $PARTITION $JOB $NCPU"
            $COMMAND > job_id_running
            # echo "Running \`$COMMAND\`"
            RUNNING_JOBS_COUNT=$(squeue -u $USER -p $PARTITION | wc -l)
            if (( $RUNNING_JOBS_COUNT >= 1000 ))
            then 
                break
            fi
        done
    done
    sleep 10m
    echo "There are currently $RUNNING_JOBS_COUNT jobs submitted to the $PARTITION cluster partition and $FINISHED_REGRESSION_COUNT jobs completed overall for $PHENOTYPE"
done
echo "Finished with regressions"


