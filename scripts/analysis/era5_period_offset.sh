#!/bin/bash
# How much of the "-2.2 K Siberian JJA cold bias" is just the reference period?
#
# eval_round10_A.py scores T2m against /work/ab0246/a270092/obs/era5/netcdf/T2M.nc,
# which is ERA5 **1990-2014**. The AMIP runs are 1870s observed SST with 1850 GHG.
# That is ~130 years of greenhouse warming folded into a number we have been
# treating as model error and tuning against.
#
# This extracts ERA5 Siberian JJA T2m per year from the DKRZ pool, 1940-2014, so the
# period offset can be measured instead of assumed. The box matches
# eval_round10_A.py (55-75N, 60-180E). NO land mask is applied -- the box is
# overwhelmingly land, and for a *difference between periods* the constant ocean
# fraction cancels. Absolute values here are therefore not directly comparable to
# the masked model numbers; the period differences are.
#
# ERA5 is on a reduced Gaussian N320 grid, which cdo cannot sellonlatbox directly,
# hence -setgridtype,regular. Operators chain right-to-left.
#
# Output: scripts/analysis/era5_siberia_jja.txt  (year, JJA box-mean K)
set -u
module load cdo 2>/dev/null
SRC=/pool/data/ERA5/E5/sf/an/1M/167
OUT="$(dirname "$0")/era5_siberia_jja.txt"
: > "$OUT"
for y in $(seq 1940 2014); do
    f="$SRC/E5sf00_1M_${y}_167.grb"
    [ -f "$f" ] || { echo "$y MISSING" >> "$OUT"; continue; }
    v=$(cdo -s -outputf,%10.5f,1 -fldmean -sellonlatbox,60,180,55,75 \
            -timmean -selmon,6,7,8 -setgridtype,regular "$f" 2>/dev/null | tr -d ' ')
    echo "$y ${v:-NaN}" >> "$OUT"
done
echo "wrote $OUT"
