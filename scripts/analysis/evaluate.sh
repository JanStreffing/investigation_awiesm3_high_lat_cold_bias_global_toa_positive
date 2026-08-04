#!/bin/bash
# THE standard evaluation. Run this, and only this, whenever a run finishes.
#
# Why a wrapper exists: the protocol was five separate scripts that had to be
# remembered and run in order, and skipping one has repeatedly produced wrong
# conclusions in this campaign -- a lever promoted on a JJA mean that was noise
# (B8), a winter penalty unnoticed for eleven rounds (B5), and a global TOA
# decomposition that was only ever done for the control. Everything that is
# supposed to be checked is now checked by default, because the way to make a
# guardrail effective is to make skipping it require effort.
#
#   usage:  ./evaluate.sh              full protocol
#           ./evaluate.sh --quick      main table + significance only
#           ./evaluate.sh --obs        also the satellite/observational checks
#
# ADDING A RUN: edit scripts/analysis/runs.py only. Every script below imports
# RUNS from there. Runs with incomplete output are skipped with a warning, so
# adding an in-flight experiment is harmless.
set -uo pipefail
cd "$(dirname "$0")"

QUICK=0; OBS=0
for a in "$@"; do
  case "$a" in
    --quick) QUICK=1 ;;
    --obs)   OBS=1 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option: $a"; exit 2 ;;
  esac
done

run () {   # run <label> <script...>
  local label="$1"; shift
  echo
  echo "################################################################################"
  echo "## $label"
  echo "################################################################################"
  python3 "$@" || echo "  !! $label FAILED (exit $?) -- continuing"
}

# 1. the main table. Includes the global TOA decomposition, which is the energy
#    half of the campaign's two targets and used to be reported as a single
#    number with no breakdown.
run "MAIN TABLE  (targets, guardrails, global TOA decomposition, RMSE)" eval_round10_A.py

# 2. significance. NEVER quote a boreal delta without the t from here.
run "NOISE FLOOR  (Siberian JJA T2m, run x year ANOVA)"                 noise_floor.py

# 3. seasons, each against ITS OWN threshold. Winter noise is 2.4x summer's.
run "SEASONAL     (per-season deltas and per-season thresholds)"        seasonal_by_run.py

if [ "$QUICK" -eq 0 ]; then
  run "RMSE SIGNIFICANCE  (per-year spatial SW RMSE, deep-water boxes)" rmse_significance.py
  run "MONTHLY LEVERS     (month-by-month T2m, SW, albedo)"             monthly_lever_check.py
fi

if [ "$OBS" -eq 1 ]; then
  # Observational checks. Slower, and they need the ERA5/satellite files from
  # albedo_decompose_prep.sh, so they are opt-in rather than default.
  run "ALBEDO DECOMPOSITION  (vs ERA5, with CERES cross-check)"         albedo_decompose.py
  run "SNOW BUDGET           (snowfall vs loss, vs ERA5)"               snow_budget.py
  run "SNOW COVER vs SATELLITE  (Rutgers/IMS melt timing)"              snowcover_vs_satellite.py
  run "BIAS BY SURFACE TYPE     (is the cold bias boreal or global?)"   bias_by_tile.py
  run "VERTICAL COLUMN          (soil L4 -> 100 hPa)"                   vertical_bias_column.py
  run "TROPOSPHERIC SECTION     (lat x pressure, and clear-sky SW)"     tropo_bias_section.py
  run "LAND ALBEDO SPLIT        (snow tile vs snow-free surface)"       land_albedo_snow_split.py
fi

cat <<'EOF'

################################################################################
## READING THE OUTPUT -- the traps that have cost this campaign real time
################################################################################
  * Quote t, or "within noise". A 4-year ranking was entirely noise once.
  * Use each season's OWN threshold: DJF +-0.588 K, MAM +-0.386, JJA +-0.242,
    SON +-0.431. Judging DJF against the JJA number overstates it by 2.4x.
  * Check tropics net TOA on anything touching convection or cloud.
  * Report absolute values, not only deltas.
  * The PI energy target is net TOA ~ 0, NOT the CERES column. CERES is
    present-day (2005-2015) and carries the real warming imbalance, so it reads
    ~+1.1 W/m2 globally. Do not tune a pre-industrial control toward it.
  * Never predict superposition -- measure it. It has been wrong in SIGN twice.
EOF
