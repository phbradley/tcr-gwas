#!/bin/zsh

type=$1
mkdir $3 
mkdir $3/${type}

for number in {01..$2..1}; do
  mkdir $3/${type}/${type}_${number}
  cd   $3/${type}/${type}_${number}
  echo "#!/bin/bash
  #SBATCH --nodes=1
  #SBATCH --mail-user=magruss@uw.edu
  #SBATCH --mail-type=END,FAIL
  source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
  conda activate r
  set -eu
  Rscript /home/mrussel2/tcr-gwas/trimming_regression/scripts/run_cluster_second_pass.R ${number} ${type} $(echo '$1') 100 PCA
  echo "done" > run.sentinel
  rm slurm*.out" > run_script_on_cluster.sh
  chmod +x run_script_on_cluster.sh
  chmod +x $PWD/run_script_on_cluster.sh
  cd /home/mrussel2/tcr-gwas/trimming_regression/send_regressions_to_cluster
done
