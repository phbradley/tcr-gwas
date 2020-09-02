#!/bin/zsh

for type in "v_trim" "j_trim" "d0_trim" "d1_trim"; do for number in {01..35500000..5000}; do
  mkdir  cluster_job_directories/${type}_${number}
  cd   cluster_job_directories/${type}_${number}
  echo "#!/bin/zsh
	#SBATCH --nodes=1
	#SBATCH --cpus-per-task=4
	#SBATCH --time=1-0
	conda activate r
	cd /home/mrussel2/tcr-gwas/trimming_regression
	Rscript run_cluster.R ${number} ${type}" > run_script_on_cluster.sh
  chmod +x run_script_on_cluster.sh
  cd /home/mrussel2/tcr-gwas/trimming_regression

done
done
