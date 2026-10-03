#!/bin/bash
#SBATCH --job-name=ihf_eval
#SBATCH --account=clidyn.clidyn
#SBATCH --partition=smp
#SBATCH --qos=30min
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --output=logs/ihf_eval_%j.log
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
module load analysis-toolbox/python-04.2026 >/dev/null 2>&1
for e in ihf0 ihf1; do
  python scripts/analysis/climate_annual_series.py $e /albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4/PICAL_crunveg_$e 2120 2129 data/clim/$e.csv > /dev/null &
  ROOT=/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4 ARM=PICAL_crunveg_$e Y0=2121 Y1=2129 python scripts/analysis/ice_albedo_requirement.py > data/clim/ice_albedo_req_$e.txt 2>&1 &
  ROOT=/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4 ARM=PICAL_crunveg_$e Y0=2121 Y1=2129 python scripts/analysis/ice_ocean_heat_by_hemisphere.py > data/clim/ice_ocean_heat_$e.txt 2>&1 &
done
wait
