#!/bin/bash
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL
source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate tcr-gwas
set -eu

Rscript /home/mrussel2/tcr-gwas/gwas_regressions/analysis_scripts/test_d0_tcrb.R $1 $2 $3
