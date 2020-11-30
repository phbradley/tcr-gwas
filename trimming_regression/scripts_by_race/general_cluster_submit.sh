#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL
source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate r
set -eu

Rscript /home/mrussel2/tcr-gwas/trimming_regression/scripts/$1
