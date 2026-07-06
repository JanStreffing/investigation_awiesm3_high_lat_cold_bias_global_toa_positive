#!/bin/bash
#SBATCH --job-name=lpjgdiag
#SBATCH --account=bb1469
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/work/bb1469/a270092/eval/lpjg_txt/lpjgdiag_%j.log
set -e
R=/work/bb1469/a270092/runtime/lpjg-spinup-develop/lpj_spinup_HR_2000y/outdata/lpj_guess
OUT=/work/bb1469/a270092/eval/lpjg_txt; mkdir -p $OUT
for v in fpc cmass dens anpp annual_burned_area; do
  first=$(ls $R/run*/$v.out 2>/dev/null | head -1)
  [ -z "$first" ] && { echo "$v: no files"; continue; }
  my=$(tail -1 "$first" | awk '{print $3}')
  head -1 "$first" > $OUT/${v}_sib.txt
  # Siberian box (Lat 55-70, Lon 60-140), last year only
  awk -v my=$my 'FNR>1 && $3==my && $2>=55 && $2<=70 && $1>=60 && $1<=140' $R/run*/$v.out >> $OUT/${v}_sib.txt
  echo "$v: year=$my rows=$(($(wc -l < $OUT/${v}_sib.txt)-1))"
done
echo "=== extraction done, summarizing ==="
source ~/loadconda.sh; conda activate reval
python3 /work/bb1469/a270092/eval/scripts/lpjg_grass_summary.py
