#!/bin/bash
#SBATCH --job-name=lpjgforc
#SBATCH --account=bb1469
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/work/bb1469/a270092/eval/lpjg_forcing_%j.log
set -e
source ~/loadconda.sh
conda activate reval
IN=/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_forcing/outdata/oifs
OUT=/work/ab0246/a270092/input/lpj-guess/oifs_forcing/AMIP_noLPJG_1d_1870-1879_TCO95_PI.nc
SC=/work/ab0246/a270092/input/lpj-guess/oifs_forcing/create_lpjg_forcing.py
echo "=== building $OUT ==="
python3 "$SC" "$IN" 1870 1879 "$OUT"
echo "=== result ==="; ls -la "$OUT"
source ~/loadconda.sh; conda activate reval
python3 -c "import xarray as xr;d=xr.open_dataset('$OUT');print('vars:',list(d.data_vars));print('dims:',dict(d.sizes))"
