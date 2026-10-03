#!/bin/bash
# Mirror the levante output needed to compare the albedo PICAL legs (2100+, merged
# LPJ-GUESS canopy work) against what came before them:
#   PICAL_ccnice 1940-2099      the same line before the LPJ-GUESS change
#   080a, 11E                   earlier coupled runs started from the CRUNCEP land state
# Only small files: LPJ-GUESS annual per-PFT tables, remapped monthly OIFS surface fields,
# FESOM's own hemispheric ice integrals, and monthly ice concentration/thickness for the
# last two decades. Layout mirrors <exp>/outdata/<model>/ so the campaign scripts work on
# it unchanged.
set -u
DST=/albedo/work/projects/p_awiesm3_cmip7/jstreffi/levante_mirror
SSH="ssh -o BatchMode=yes"
mkdir -p $DST

pull() {  # pull <remote exp root> <local exp name> <filter file>
    mkdir -p $DST/$2
    rsync -a --prune-empty-dirs -e "$SSH" --include-from="$3" \
        "levante:$1/outdata/" "$DST/$2/outdata/"
    echo "$(date +%F_%T) done $2 rc=$?"
}

F=$(mktemp)
cat > $F <<'EOF'
+ */
+ lpj_guess/*/run1/fpc.out
+ lpj_guess/*/run1/lai.out
+ lpj_guess/*/run1/anpp.out
+ lpj_guess/*/run1/dens.out
+ lpj_guess/*/run1/cmass.out
+ lpj_guess/*/run1/est_limits.out
+ oifs/atm_remapped_1m_2t_*.nc
+ oifs/atm_remapped_1m_cvh_*.nc
+ oifs/atm_remapped_1m_cvl_*.nc
+ oifs/atm_remapped_1m_lai_hv_*.nc
+ oifs/atm_remapped_1m_lai_lv_*.nc
+ oifs/atm_remapped_1m_ci_*.nc
+ oifs/atm_remapped_1m_fal_*.nc
+ oifs/atm_remapped_1m_sd_*.nc
+ oifs/atm_remapped_1m_skt_*.nc
+ oifs/atm_remapped_1m_tsr_*.nc
+ oifs/atm_remapped_1m_ttr_*.nc
+ oifs/atm_remapped_1m_ssr_*.nc
+ oifs/atm_remapped_1m_str_*.nc
+ oifs/atm_remapped_1m_slhf_*.nc
+ oifs/atm_remapped_1m_sshf_*.nc
+ oifs/atm_remapped_1m_sf_*.nc
+ oifs/atm_remapped_1m_tcc_*.nc
+ fesom/si*.fesom.*.nc
- fesom/*.gr.*
+ fesom/a_ice.fesom.20[89]?.nc
+ fesom/m_ice.fesom.20[89]?.nc
+ fesom/m_snow.fesom.20[89]?.nc
- *
EOF

pull /work/bb1469/a270092/runtime/awiesm3-v3.4/PICAL_ccnice PICAL_ccnice $F
pull /work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_080a_lpjguess_Baseline_coupled_fromCRUNCEP 080a $F
E11=$($SSH levante 'ls -d /work/bb1469/a270092/runtime/*/Tuning_test_11E_swemin15_K1 2>/dev/null | head -1' 2>/dev/null)
echo "11E root: $E11"
[ -n "$E11" ] && pull $E11 11E $F
rm -f $F
du -sh $DST/*
