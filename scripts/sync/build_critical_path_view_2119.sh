#!/bin/bash
# Outdata tree that reads as the critical path of the tuning campaign up to the last run before the
# hemispheric GM change, 1850-2119, for reval's spin-up plots (release_evaluation_tool2 branch
# feat/critical-path-spinup; config awiesm3_v35_critical_path_1850-2119_albedo.py). Symlinks only.
# Same construction as build_critical_path_view.sh, stopping at PICAL_crunveg (2100-2119).
set -e
P=/albedo/work/projects/p_awiesm3_cmip7/jstreffi
V=$P/reval/views/critical_path_1850-2119
C=$P/reval/spinup_cache
R=$P/runtime/awiesm3-v3.4
rm -rf $V; mkdir -p $V/outdata/oifs $V/outdata/fesom $V/restart/oasis3mct
ln -s $R/PICAL_crunveg/restart/oasis3mct/areas.nc $V/restart/oasis3mct/areas.nc
ln -s $R/PICAL_crunveg/config $V/config
ln -s $R/PICAL_crunveg/run_21100101-21191231 $V/run_21100101-21191231
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
echo "oifs links: $(ls $V/outdata/oifs | wc -l)  (expect $((270*8)))"
echo "a_ice years: $(ls $V/outdata/fesom | grep -c '^a_ice\.fesom')  temp years: $(ls $V/outdata/fesom | grep -c '^temp\.fesom')"
K=$C/critical_path_2119
rm -rf $K; mkdir -p $K
for y in $(seq 1850 2119); do for f in $C/critical_path/*_${y}.npy; do ln -s $f $K/; done; done
echo "cache files: $(ls $K | wc -l)  (expect $((270*2)))"
