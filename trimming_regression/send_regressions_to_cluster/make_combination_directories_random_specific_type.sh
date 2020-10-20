#!/bin/zsh

type=$1
mkdir random_effects_cluster_job_directory_100_boots_subset
mkdir random_effects_cluster_job_directory_100_boots_subset/${type}

for number in {01..$2..1}; do
  mkdir random_effects_cluster_job_directory_100_boots_subset/${type}/${type}_${number}
  cd   random_effects_cluster_job_directory_100_boots_subset/${type}/${type}_${number}
  echo "#!/bin/bash
  #SBATCH --nodes=1
  #SBATCH --mail-user=magruss@uw.edu
  #SBATCH --mail-type=END,FAIL
  source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
  conda activate r
  set -eu
  Rscript /home/mrussel2/tcr-gwas/trimming_regression/scripts/run_cluster_random_by_portion_of_file.R ${number} ${type} $(echo '$1') 100
  echo "done" > run.sentinel
  rm slurm*.out" > run_script_on_cluster.sh
  chmod +x run_script_on_cluster.sh
  chmod +x $PWD/run_script_on_cluster.sh
  cd /home/mrussel2/tcr-gwas/trimming_regression/send_regressions_to_cluster
done
