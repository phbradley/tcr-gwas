#!/bin/sh

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate tcr-gwas
set -eu

PHENOTYPE=$1
NCPU=$2

Rscript /home/mrussel2/tcr-gwas/gwas_regressions/analysis_scripts/bootstrap_analysis.R $1 $2

