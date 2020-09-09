#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL
source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu
cpulimit -l 400 -i Rscript /home/mrussel2/tcr-gwas/trimming_regression/run_cluster.R ${number} ${type}
echo "done" > run.sentinel
rm slurm*.out
