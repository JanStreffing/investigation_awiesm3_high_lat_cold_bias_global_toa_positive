#!/bin/bash
# One outdata tree that reads as the critical path of the tuning campaign, 1850-2139, for reval's
# spin-up plots (release_evaluation_tool2 branch feat/critical-path-spinup, config
# awiesm3_v35_critical_path_1850-2139_albedo.py). Symlinks only.
#   1850-1919 PICAL, 1920-1939 PICAL_momixoff   levante_mirror/PICAL_momixoff_1850-1939 (momixoff links PICAL's years)
#   1940-2099 PICAL_ccnice                      levante_mirror/PICAL_ccnice
#   2100-2119 PICAL_crunveg, 2120-2139 PICAL_crunveg_gmhemi1800   on albedo
# Before 2100 only the monthly OpenIFS surface fields and a_ice exist here; the 3D ocean diagnostics
# for those years come from reval's spin-up cache, computed on levante.
# tsrc/ttrc are deliberately left out: the ccnice mirror does not have them, and part2 needs a
# variable for every year or not at all.
set -e
P=/albedo/work/projects/p_awiesm3_cmip7/jstreffi
V=$P/reval/views/critical_path_1850-2139
R=$P/runtime/awiesm3-v3.4
rm -rf $V; mkdir -p $V/outdata/oifs $V/outdata/fesom $V/restart/oasis3mct
ln -s $R/PICAL_crunveg_gmhemi1800/restart/oasis3mct/areas.nc $V/restart/oasis3mct/areas.nc
# reval's mask plot and cavity check look next to outdata/ for run_*/work/masks.nc and config/fesom/namelist.config
ln -s $R/PICAL_crunveg_gmhemi1800/config $V/config
ln -s $R/PICAL_crunveg_gmhemi1800/run_21200101-21391231 $V/run_21200101-21391231
link_years() {   # <source outdata> <first> <last> <mode: light|full>
  local src=$1 y
  for y in $(seq $2 $3); do
    for v in 2t ssr str tsr ttr sf slhf sshf; do
      f=$src/oifs/atm_remapped_1m_${v}_${y}-${y}.nc; [ -e $f ] && ln -s $f $V/outdata/oifs/ || echo "MISSING $f"
    done
    if [ "$4" = full ]; then
      for f in $src/fesom/*.fesom.$y.nc; do case $f in *.gr.*) ;; *) ln -s $f $V/outdata/fesom/ ;; esac; done
    else
      f=$src/fesom/a_ice.fesom.$y.nc; [ -e $f ] && ln -s $f $V/outdata/fesom/ || echo "MISSING $f"
    fi
  done
}
link_years $P/levante_mirror/PICAL_momixoff_1850-1939/outdata 1850 1939 light
link_years $P/levante_mirror/PICAL_ccnice/outdata            1940 2099 light
link_years $R/PICAL_crunveg/outdata                           2100 2119 full
link_years $R/PICAL_crunveg_gmhemi1800/outdata                2120 2139 full
echo "oifs links: $(ls $V/outdata/oifs | wc -l)  (expect $((290*8)))"
echo "a_ice years: $(ls $V/outdata/fesom | grep -c '^a_ice\.fesom')  temp years: $(ls $V/outdata/fesom | grep -c '^temp\.fesom')  w years: $(ls $V/outdata/fesom | grep -c '^w\.fesom')"
