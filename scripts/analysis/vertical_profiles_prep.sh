#!/bin/bash
# Build matched model/ERA5 pressure-level climatologies for vertical profile comparison.
#
# WHY. Everything the campaign has measured so far is at the surface or TOA. That tells us a
# bias exists but not WHERE in the column it originates -- e.g. whether the boreal JJA cold
# bias is confined to the boundary layer (a surface-exchange problem) or extends through the
# troposphere (a circulation or radiation problem), and whether the Southern Ocean cloud
# deficit sits in low cloud (consistent with the D2b pressure-gate result) or aloft.
#
# TWO COMPARISONS:
#   PRIMARY   amip_presentday (1990-2014) vs ERA5 (1990-2014)  -- period-clean, no epoch offset
#   SECONDARY amip_pi_base    (1872-1915) vs ERA5 (1940-1969)  -- closest ERA5 epoch to the runs
# The primary is the trustworthy one; the secondary carries the residual epoch offset that
# amip_presentday measured as +0.42 K on Siberian JJA T2m.
#
# DATA. Model writes MONTHLY pressure-level output (36 MB/yr/var) as well as 6-hourly (4.6 GB),
# so this uses the monthly files. Variables available: t, q, r, u, v, w, vo, z.
# NOTE: the model does NOT write cloud variables on pressure levels (no cc/clwc/ciwc), so
# vertical CLOUD structure cannot be compared directly -- RH is the proxy.
#
# ERA5 comes from the DKRZ pool (no download needed): /pool/data/ERA5/E5/pl/an/1M/<param>/
#   130 = temperature, 157 = relative humidity, 133 = specific humidity
# It is on a reduced Gaussian N320 grid with 37 levels; all 19 model levels exist exactly in
# that set, so -sellevel matches with NO interpolation. -setgridtype,regular is required
# before any regridding, and -f nc is required or cdo writes GRIB.
#
# Output: /tmp/vprof/{era5,model}_<tag>_<var>_<season>.nc  on the model grid, 19 levels.
set -u
module load cdo 2>/dev/null
# NOT /tmp: it is cleaned periodically and this cache was lost twice mid-analysis,
# once taking an hour of ERA5 pool reads with it.  /work persists.
W=${VPROF_DIR:-/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/data/vprof}; mkdir -p $W
RT=/work/bb1469/a270092/runtime/oifsamip-cy48
E5=/pool/data/ERA5/E5/pl/an/1M
LEV=100000,92500,85000,70000,60000,50000,40000,30000,25000,20000,15000,10000,7000,5000,3000,2000,1000,500,100
declare -A PARAM=( [t]=130 [r]=157 [q]=133 )
declare -A SEASON=( [JJA]=6,7,8 [DJF]=12,1,2 [ANN]=1,2,3,4,5,6,7,8,9,10,11,12 )

# grid description from any model file
cdo -s griddes $(ls $RT/amip_presentday/outdata/oifs/atm_remapped_1m_pl_t_1m_pl_*.nc | head -1) > $W/modelgrid.txt

do_era5 () {  # tag y0 y1
  local tag=$1 y0=$2 y1=$3
  for v in t r q; do
    for s in JJA DJF ANN; do
      local out=$W/era5_${tag}_${v}_${s}.nc
      [ -f "$out" ] && continue
      local yrs=""
      for y in $(seq $y0 $y1); do
        f=$E5/${PARAM[$v]}/E5pl00_1M_${y}_${PARAM[$v]}.grb
        [ -f "$f" ] || continue
        cdo -s -O -f nc -remapbil,$W/modelgrid.txt -sellevel,$LEV -setgridtype,regular \
            -timmean -selmon,${SEASON[$s]} "$f" $W/.e5_${v}_${s}_${y}.nc 2>/dev/null && yrs="$yrs $W/.e5_${v}_${s}_${y}.nc"
      done
      [ -n "$yrs" ] && cdo -s -O ensmean $yrs "$out" && rm -f $yrs
      echo "  era5 $tag $v $s -> $(basename $out)"
    done
  done
}

do_model () {  # tag exp y0 y1
  local tag=$1 exp=$2 y0=$3 y1=$4
  for v in t r q; do
    for s in JJA DJF ANN; do
      local out=$W/model_${tag}_${v}_${s}.nc
      [ -f "$out" ] && continue
      local yrs=""
      for y in $(seq $y0 $y1); do
        f=$RT/$exp/outdata/oifs/atm_remapped_1m_pl_${v}_1m_pl_${y}-${y}.nc
        [ -f "$f" ] || continue
        cdo -s -O -timmean -selmon,${SEASON[$s]} "$f" $W/.m_${v}_${s}_${y}.nc 2>/dev/null && yrs="$yrs $W/.m_${v}_${s}_${y}.nc"
      done
      [ -n "$yrs" ] && cdo -s -O ensmean $yrs "$out" && rm -f $yrs
      echo "  model $tag $v $s -> $(basename $out)"
    done
  done
}

echo "PRIMARY: present-day, period-clean"
do_era5  pd 1990 2014
do_model pd amip_presentday 1990 2014
echo "SECONDARY: PI epoch"
do_era5  pi 1940 1969
do_model pi amip_pi_base 1872 1915
echo "done"
