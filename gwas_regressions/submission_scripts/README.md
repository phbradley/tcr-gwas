# Submission scripts

This directory contains files for running the analyses genome-wide for a variety of phenotypes GENOME-WIDE by splitting up the regressions into cluster jobs.
The files are as follows:

* [Create a directory for each cluster job](make_job_directories.sh) which contains the appropriate scripts for each job within each job directory
* [Create job submission file](make_submission_file.sh) for a specified phenotype class. We also have similar scripts for creating job sumbmission files for the [N-insertion class](make_submission_file_inserts.sh) and for a [single phenotype](make_submission_file_single.sh).
* [Submit jobs to the cluster continuously](submit_cluster_jobs_continuously.sh) to a specified cluster partition and node count, for a specified phenotype. 
* All of the above files are sourced by the top level [genome-wide analysis submission script](run_regressions_genome_wide.sh)
* A script to run a series of random [bootstrap analyses](bootstrap_analysis.sh) (for downstream lambda calculations)

Phenotype specific job submission functions and variables are located in the following [directory](phenotype_functions).


