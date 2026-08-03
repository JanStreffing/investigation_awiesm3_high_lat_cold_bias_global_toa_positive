#!/bin/bash
# ERA5 snow/albedo fields for the June albedo decomposition, 1990-2014.
#
# Pulls the four fields that let the surface albedo be split into cover, snow
# brightness and snow-free brightness:
#   243 fal  forecast (surface) albedo
#   032 asn  snow albedo
#   141 sd   snow depth (water equivalent, m)
#   033 rsn  snow density (kg/m3)
#
# ERA5 assimilates snow-cover observations (IMS), so its snow EXTENT is
# observationally constrained -- which is what makes it a usable reference for
# "does the model keep snow too long in June", even though its albedo scheme is
# a relative of HTESSEL and so is NOT independent for the brightness terms.
#
# Traps handled here (see .claude/skills/amip-eval/SKILL.md):
#   * ERA5 is on a reduced Gaussian grid -> -setgridtype,regular BEFORE anything
#   * cdo writes GRIB unless told otherwise -> -f nc
#   * cdo operator chains read RIGHT TO LEFT
set -euo pipefail

POOL=/pool/data/ERA5/E5/sf/an/1M
OUT=${1:-/work/ab0246/a270092/obs/era5/snow}
Y0=1990; Y1=2014

mkdir -p "$OUT"
module load cdo 2>/dev/null || true

# analysis stream: 243 fal, 032 asn, 141 sd, 033 rsn, 238 tsn (snow temperature)
# forecast stream (accumulated fluxes, under sf/fc/1M):
#   144 sf   snowfall      045 smlt snowmelt
#   146 sshf sensible HF   147 slhf latent HF   (positive DOWNWARD, into surface)
#   176 ssr  net solar     177 str  net thermal
#   175 strd down thermal  169 ssrd down solar
for p in 243 032 141 033 238; do
  tgt="$OUT/era5_${p}_clim_${Y0}-${Y1}.nc"
  if [ -f "$tgt" ]; then echo "have $(basename "$tgt")"; continue; fi
  echo "building $(basename "$tgt") ..."
  files=()
  for y in $(seq $Y0 $Y1); do
    f="$POOL/$p/E5sf00_1M_${y}_${p}.grb"
    [ -f "$f" ] && files+=("$f")
  done
  if [ ${#files[@]} -eq 0 ]; then echo "  !! no input for $p"; continue; fi
  # monthly climatology over the period, on a regular grid, as netCDF
  cdo -s -f nc -ymonmean -setgridtype,regular -cat "${files[@]}" "$tgt"
done

echo
echo "done. contents:"
for f in "$OUT"/era5_*_clim_${Y0}-${Y1}.nc; do
  [ -f "$f" ] && printf "  %-40s %s\n" "$(basename "$f")" "$(cdo -s showname "$f" 2>/dev/null | tr -s ' ')"
done

# --- snow BUDGET fluxes: these are forecast (accumulated) fields, different stream
POOLF=/pool/data/ERA5/E5/sf/fc/1M
for p in 144 045 146 147 176 177 175 169; do
  tgt="$OUT/era5_${p}_clim_${Y0}-${Y1}.nc"
  if [ -f "$tgt" ]; then echo "have $(basename "$tgt")"; continue; fi
  echo "building $(basename "$tgt") ..."
  files=()
  for y in $(seq $Y0 $Y1); do
    f="$POOLF/$p/E5sf12_1M_${y}_${p}.grb"
    [ -f "$f" ] && files+=("$f")
  done
  if [ ${#files[@]} -eq 0 ]; then echo "  !! no input for $p"; continue; fi
  cdo -s -f nc -ymonmean -setgridtype,regular -cat "${files[@]}" "$tgt"
done
