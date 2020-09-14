#!/bin/zsh

for type in "v_trim" "j_trim" "d0_trim" "d1_trim"; do for number in {01..10000000..1000}; do
  cd   random_effects_cluster_job_directory/${type}/${type}_${number}
  echo "#!/bin/bash
  #SBATCH --nodes=1
  #SBATCH --mail-user=magruss@uw.edu
  #SBATCH --mail-type=END,FAIL
  source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
  conda activate r
  set -eu
  cpulimit -l 400 -i Rscript /home/mrussel2/tcr-gwas/trimming_regression/run_cluster_random.R ${number} ${type} $(echo '$1')
  echo "done" > run.sentinel
  rm slurm*.out" > run_script_on_cluster.sh
  chmod +x run_script_on_cluster.sh
  chmod +x $PWD/run_script_on_cluster.sh
  cd /home/mrussel2/tcr-gwas/trimming_regression

done
done
