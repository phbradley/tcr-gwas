#!/bin/zsh
	#SBATCH --nodes=1
	#SBATCH --cpus-per-task=4
	#SBATCH --time=1-0
	conda activate r
	cd ../../
	Rscript run_cluster.R 32650001 v_trim
