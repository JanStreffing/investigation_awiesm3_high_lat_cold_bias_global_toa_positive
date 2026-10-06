#!/bin/bash
# Outdata trees that read as the critical path of the tuning campaign continued to 2209, one per
# candidate line, for reval's spin-up plots (release_evaluation_tool2 branch
# feat/critical-path-spinup; configs awiesm3_v35_critical_path_{tke,idemix}_1850-2209_albedo.py).
# Symlinks only. Same construction as build_critical_path_view.sh (1850-2139), then
#   2140-2159 PICAL_crunveg_gmhemi_tke          TKE only (IDEMIX off)
#   tke:    2160-2209 PICAL_crunveg_tke_albsn082      albsn 0.82
#   idemix: 2160-2169 PICAL_crunveg_tke_albsn082, 2170-2209 PICAL_crunveg_idemix_albsn082 (IDEMIX on again)
# The spin-up cache (AMOC, Hovmoeller) is assembled the same way from the per-run caches.
# Usage: build_critical_path_view_2209.sh tke|idemix
set -e
L=$1
P=/albedo/work/projects/p_awiesm3_cmip7/jstreffi
V=$P/reval/views/critical_path_${L}_1850-2209
C=$P/reval/spinup_cache
R=$P/runtime/awiesm3-v3.4
case $L in
  tke)    LAST=PICAL_crunveg_tke_albsn082 ;;
  idemix) LAST=PICAL_crunveg_idemix_albsn082 ;;
  *) echo "usage: $0 tke|idemix"; exit 1 ;;
esac
rm -rf $V; mkdir -p $V/outdata/oifs $V/outdata/fesom $V/restart/oasis3mct
ln -s $R/$LAST/restart/oasis3mct/areas.nc $V/restart/oasis3mct/areas.nc
# reval's mask plot and cavity check look next to outdata/ for run_*/work/masks.nc and config/fesom/namelist.config
ln -s $R/$LAST/config $V/config
ln -s $R/$LAST/run_21900101-22091231 $V/run_21900101-22091231
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
link_years $R/PICAL_crunveg/outdata                           2100 2119 light
link_years $R/PICAL_crunveg_gmhemi1800/outdata                2120 2139 light
link_years $R/PICAL_crunveg_gmhemi_tke/outdata                2140 2159 light
if [ $L = tke ]; then
  link_years $R/PICAL_crunveg_tke_albsn082/outdata            2160 2209 full
else
  link_years $R/PICAL_crunveg_tke_albsn082/outdata            2160 2169 light
  link_years $R/PICAL_crunveg_idemix_albsn082/outdata         2170 2209 full
fi
echo "oifs links: $(ls $V/outdata/oifs | wc -l)  (expect $((360*8)))"
echo "a_ice years: $(ls $V/outdata/fesom | grep -c '^a_ice\.fesom')  temp years: $(ls $V/outdata/fesom | grep -c '^temp\.fesom')"

K=$C/critical_path_${L}
rm -rf $K; mkdir -p $K
cache_years() {   # <cache dir> <first> <last>
  local y f
  for y in $(seq $2 $3); do for f in $1/*_${y}.npy; do ln -s $f $K/; done; done
}
cache_years $C/critical_path 1850 2139
cache_years $C/gmhemi_tke    2140 2159
if [ $L = tke ]; then cache_years $C/tke_albsn082 2160 2209
else cache_years $C/tke_albsn082 2160 2169; cache_years $C/idemix_albsn082 2170 2209; fi
echo "cache files: $(ls $K | wc -l)  (expect $((360*2)))"

# LPJ-GUESS view of the last run: reval expects <range>/run1/
G=$P/reval/views/$LAST/lpj_guess
for d in $R/$LAST/outdata/lpj_guess/*/; do r=$(basename $d); mkdir -p $G/$r; ln -sfn ${d%/} $G/$r/run1; done
ls $G
