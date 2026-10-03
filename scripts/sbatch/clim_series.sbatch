#!/bin/bash
#SBATCH --job-name=clim_series
#SBATCH --account=clidyn.clidyn
#SBATCH --partition=smp
#SBATCH --qos=30min
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/clim_series_%j.log
# usage: sbatch clim_series.sbatch tag=root:y0:y1 [...]
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
module load analysis-toolbox/python-04.2026 >/dev/null 2>&1
for a in "$@"; do
  tag=${a%%=*}; rest=${a#*=}; root=${rest%%:*}; yy=${rest#*:}; y0=${yy%%:*}; y1=${yy#*:}
  python scripts/analysis/climate_annual_series.py $tag $root $y0 $y1 data/clim/${tag}.csv &
done
wait
