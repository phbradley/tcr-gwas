#!/bin/bash
  #SBATCH --nodes=1
  #SBATCH --mail-user=magruss@uw.edu
  #SBATCH --mail-type=END,FAIL
  source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
  conda activate r
  set -eu
  Rscript /home/mrussel2/tcr-gwas/trimming_regression/scripts/run_cluster_random_by_portion_of_file.R 178 v_trim $1 100
  echo done > run.sentinel
  rm slurm*.out
