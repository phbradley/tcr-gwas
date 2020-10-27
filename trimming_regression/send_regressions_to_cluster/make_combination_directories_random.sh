#!/bin/zsh

mkdir $1
mkdir $1/d1_trim
mkdir $1/d0_trim
mkdir $1/v_trim
mkdir $1/j_trim
mkdir $1/vd_insert
mkdir $1/vj_insert
mkdir $1/dj_insert

for type in "v_trim" "d1_trim" "d0_trim" "j_trim" "vd_insert" "vj_insert" "dj_insert"; do for number in {01..35481497..1000}; do
  mkdir $1/${type}/${type}_${number}
  cd  $1/${type}/${type}_${number}
  echo "#!/bin/bash
  #SBATCH --nodes=1
  #SBATCH --mail-user=magruss@uw.edu
  #SBATCH --mail-type=END,FAIL
  source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
  conda activate r
  set -eu
  Rscript /home/mrussel2/tcr-gwas/trimming_regression/scripts/run_cluster_random.R ${number} ${type} $(echo '$1') 0
  echo "done" > run.sentinel
  rm slurm*.out" > run_script_on_cluster.sh
  chmod +x run_script_on_cluster.sh
  chmod +x $PWD/run_script_on_cluster.sh
  cd /home/mrussel2/tcr-gwas/trimming_regression/send_regressions_to_cluster
done
done
