#!/bin/bash
#SBATCH --job-name=lpjg2nc
#SBATCH --account=bb1469
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=/work/bb1469/a270092/eval/lpjg_nc/lpjg2nc_%j.log
set -e
source ~/loadconda.sh
conda activate reval
BRIDGE=/work/bb1469/a270092/eval/lpjg_bridge
OUT=/work/bb1469/a270092/eval/lpjg_nc
SC=/work/ab0246/a270092/software/lpjg2nc2/lpjg2nc.py
# key discriminators first: dens (establishment vs mortality), burned area (fire), then cmass/anpp/fpc
for f in dens.out annual_burned_area.out cmass.out anpp.out fpc.out; do
  echo "=== converting $f ($(date +%T)) ==="
  python $SC -p "$BRIDGE" -f "$f" -o "$OUT" -j 8 --inner-jobs 16 && echo "OK: $f" || echo "FAILED: $f"
done
echo "=== DONE $(date +%T) ==="; ls -la $OUT/*.nc
