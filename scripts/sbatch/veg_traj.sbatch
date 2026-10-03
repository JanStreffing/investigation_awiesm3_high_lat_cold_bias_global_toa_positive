#!/bin/bash
#SBATCH --job-name=veg_traj
#SBATCH --account=clidyn.clidyn
#SBATCH --partition=smp
#SBATCH --qos=30min
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:30:00
#SBATCH --output=logs/veg_traj_%j.log
# usage: sbatch veg_traj.sbatch tag=root [tag=root ...]
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
module load analysis-toolbox/python-04.2026 >/dev/null 2>&1
for a in "$@"; do
  tag=${a%%=*}; root=${a#*=}
  for box in sib boreal; do
    python scripts/analysis/boreal_veg_trajectory.py $tag $root data/veg_traj/${tag}_${box}.csv $box &
  done
done
wait
