#!/bin/bash
# Does the ENERGY target carry the same reference-period flaw as the boreal one?
#
# The SO SW CRE gap (+8.08 W/m2) and the global TOA imbalance (+0.53 W/m2) are scored
# against CERES EBAF Ed4.1, whose climatology here is 07/2005-06/2015. The AMIP runs are
# 1870s observed SST with 1850 GHG. Same ~130-year mismatch that turned out to be half of
# the boreal "bias" -- but this time on the target A1b was tuned against, and never checked.
#
# CERES itself cannot answer this: it starts in 2000. ERA5 can, from 1940, using the same
# quantities:
#     SO SW CRE   = TSR - TSRC  (178 - 208) over 45-65S
#     global TOA  = TSR + TTR   (178 + 179)
#
# IMPORTANT CAVEAT, stated up front: ERA5's TOA radiative fluxes are *model-derived*
# forecast fields, not assimilated observations. Their long-term drift partly reflects
# changes in the observing system and in the forecast model, so these numbers are a
# plausibility estimate of the period offset -- far weaker evidence than the T2m offset,
# which rests on an assimilated variable. Treat as an order-of-magnitude bound.
#
# No ocean mask is applied to the 45-65S band: it is ~97% ocean, and the same band is used
# for both periods so the small land fraction cancels in the difference.
#
# Output: scripts/analysis/era5_energy_periods.txt   (year  so_swcre  global_toa)
set -u
module load cdo 2>/dev/null
P=/pool/data/ERA5/E5/sf/fc/1M
OUT="$(dirname "$0")/era5_energy_periods.txt"
: > "$OUT"
# ERA5 monthly-mean fluxes are accumulations in J/m2 over the accumulation period;
# divide by seconds to get W/m2. For 1M means of 12-hourly accumulations the step is 86400 s.
for y in $(seq 1940 2015); do
    a=$P/178/E5sf12_1M_${y}_178.grb; b=$P/208/E5sf12_1M_${y}_208.grb; c=$P/179/E5sf12_1M_${y}_179.grb
    [ -f "$a" ] && [ -f "$b" ] && [ -f "$c" ] || { echo "$y MISSING" >> "$OUT"; continue; }
    cre=$(cdo -s -outputf,%12.5f,1 -fldmean -sellonlatbox,-180,180,-65,-45 -timmean \
              -sub -setgridtype,regular "$a" -setgridtype,regular "$b" 2>/dev/null | tr -d ' ')
    toa=$(cdo -s -outputf,%12.5f,1 -fldmean -timmean \
              -add -setgridtype,regular "$a" -setgridtype,regular "$c" 2>/dev/null | tr -d ' ')
    echo "$y ${cre:-NaN} ${toa:-NaN}" >> "$OUT"
done
echo "wrote $OUT"
