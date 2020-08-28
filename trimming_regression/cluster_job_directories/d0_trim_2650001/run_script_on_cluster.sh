#!/bin/zsh
	#SBATCH --nodes=1
	#SBATCH --cpus-per-task=4
	#SBATCH --time=1-0
	conda activate r
	cd ../../
	Rscript run_cluster.R 2650001 d0_trim
