# Data paths & experiments

All paths are on DKRZ **Levante**. `a270092` and `a270270` are the two user
accounts holding the runs; storage is under `/work/bb1469` (compute project) and
`/work/ab0246` (input/software).

## Experiments

### 1. Imbalance-tuning ensemble (Part I) — 22 coupled AWI-ESM3 runs
CORE3 mesh / TCO95 atmosphere / LPJ-GUESS, 30-yr PI-control segments (model years
1350–1379), all branched from `Tuning_test_06_Baseline`.

| Run | Path |
|---|---|
| Baseline + all `06*` variants | `/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06*` |

Per-run distinguishing parameter change: report App. A (or `report/report.tex`).
Per-run `reval` diagnostics staged under `run_figures/<run>/` (radiation_budget,
t2m_vs_ERA5, cmpi).

### 2. Clean no-LPJG CORE3 reference (Part II — the WITH/WITHOUT pair)
Same CORE3 ocean + TCO95 atmosphere as Baseline but **HTESSEL + prescribed
satellite vegetation** (no interactive LPJG). Differs from Baseline in LPJG alone.

```
/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y
  outdata/oifs/atm_remapped_1m_<var>_1m_<yr>-<yr>.nc      # monthly remapped fields
```
The coupled WITH-LPJG side of the pair is `Tuning_test_06_Baseline` (above).

### 3. HTESSEL LPJG spin-up (dominant-PFT reference)
```
/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_54_100YRES_2000YSPINUP_TCO95_CORE3
```

### 4. OIFS–AMIP forcing generator (Part II decisive test — feedback-free climate)
Clean OpenIFS-AMIP TCO95 run, **no LPJG**, prescribed observed SST (~PI 1870–79),
emitting the daily fields LPJG needs as offline forcing.
```
run:      /work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_forcing/outdata/oifs/atm_1d_<yr>-<yr>.nc
forcing:  /work/ab0246/a270092/input/lpj-guess/oifs_forcing/
            AMIP_noLPJG_1d_1870-1879_TCO95_PI_fluxfix.nc   # <-- USE THIS (flux-corrected)
            AMIP_noLPJG_1d_1870-1879_TCO95_PI.nc           # uncorrected (fluxes 6x too low; see report §12.2)
```

### 5. Decisive standalone LPJG re-spin (v3) — the climate-vs-competition result
2000-yr offline LPJG spin-up on the corrected AMIP forcing, dense 11 538-cell
gridlist. **This is the run behind `plots/v3_decisive_climate_vs_veg.png`.**
```
/work/bb1469/a270092/runtime/lpjg-spinup-develop/LR_2000y_PIforcing_v3
  outdata/lpj_guess/19000101-38991231/run1/*.out          # combined per-variable (fpc, anpp, cmass, dens, lai, ...)
```
Earlier, superseded spin-ups (documented for provenance):
```
.../LR_2000y_PIforcing        # chunked, restart-OOM (abandoned)
.../LR_2000y_PIforcing_1job   # single job, bad-ins (freenyears>1000)
.../LR_1000y_spinup           # old sparse/cold-sampled LR spin-up (inconclusive)
```

### 6. Coupled CRUNCEP3-initialized branch follow-up (2026-07-24) — 50-year regional land/ocean split check

Coupled AWI-ESM3 + LPJ-GUESS run started from the CRUNCEP3-based vegetation state
and followed for 50 model years. Used to test whether the boreal-forest problem
remains purely a vegetation-competition issue once the system is allowed to
adjust in coupled mode.

```
/work/bb1469/a270270/runtime/awiesm3-v3.4/<CRUNCEP3-initialized-coupled-branch-directory>
  run_13500101-13591231/work/run*/output/fpc.out   # LPJ shards for year 1350
  run_13600101-13691231/work/run*/output/fpc.out   # LPJ shards for year 1360
  run_13900101-13991231/work/run*/output/fpc.out   # LPJ shards for year 1399
  outdata/oifs/atm_remapped_1m_2t_1350-1359.nc     # via per-year files atm_remapped_1m_2t_1350-1350.nc ... 1359-1359.nc
  outdata/oifs/atm_remapped_1m_2t_1390-1399.nc     # via per-year files atm_remapped_1m_2t_1390-1390.nc ... 1399-1399.nc
```

Actual directory name on disk (internal alias in quotes):

```
"Tuning_test_080a_lpjguess_Baseline_coupled_fromCRUNCEP"
```

Generated diagnostics provenance (original generation directory):

```
/work/ab0995/a270270/analysis/LR_offline_LPJ_GUESS_tunning_new/plots/diff_from_amip_and_era5/
  treefpc_coupled_cruncep3_init_y10_minus_y0_nh.png
  treefpc_coupled_cruncep3_init_y50_minus_y0_nh.png
  t2m_seasonal_bias_coupled_cruncep3_init_vs_cruncep3_nh.png
  t2m_seasonal_bias_coupled_cruncep3_init_1390_1399_vs_cruncep3_nh.png
  t2m_bias_vs_base06_all_tuned_coupled_nh_10yrseasonal_DJF.png
  t2m_bias_vs_base06_all_tuned_coupled_nh_10yrseasonal_MAM.png
  t2m_bias_vs_base06_all_tuned_coupled_nh_10yrseasonal_JJA.png
  t2m_bias_vs_base06_all_tuned_coupled_nh_10yrseasonal_SON.png
  net_radiation_panel_baseline_A_to_V_07A_to_C_HRgoal_LRpi_last30.png
```

The same figure set has been copied into this investigation repository at:

```
/work/ab0995/a270270/analysis/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/plots/
```

Working interpretation and proposed next runs:

```
/work/ab0995/a270270/analysis/LR_offline_LPJ_GUESS_tunning_new/next_experiment_design_coupled_cruncep3_init_followup.md
```

### 7. Four standalone forcing-sensitivity spin-ups (2026-07-29)

All four are 2000-year TCO95 offline LPJ-GUESS spin-ups. Final-year diagnostics
use model year 3900 from `outdata/lpj_guess/19000101-38991231/run1/`.

| Experiment | Run directory | Forcing source |
|---|---|---|
| AMIP | `/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_2000Y_TCO95_CORE3_AMIPforcing` | `/work/ab0995/a270270/input/cruncep_v7/AMIP_noLPJG_1d_1870-1879_TCO95_PI_fluxfix.nc` |
| CRUNCEP v3 | `/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_2000Y_TCO95_CORE3_CRUNCEPcalibrated_v3` | `/work/ab0995/a270270/input/cruncep_v7/CRUNCEP_noLPJG_1d_1901-1910_TCO95_calibrated_v3.nc` |
| CRUNCEP + CERES direct | `/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_2000Y_TCO95_CORE3_CRUNCEPandCERES` | `/work/ab0995/a270270/input/cruncep_v7/CRUNCEP_noLPJG_1d_1901-1910_TCO95_CERES_rsns_rlns.nc` |
| CRUNCEP + CERES daily variability | `/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_2000Y_TCO95_CORE3_CRUNCEPandCERES_daily_variability` | `/work/ab0995/a270270/input/cruncep_v7/CRUNCEP_noLPJG_1d_1901-1910_TCO95_CERESclim_daily_variability.nc` |

The two runtime files named
`run_19000101-38991231/work/ifs_spinup_forcing/AM04_atm_cmip6_1d_1990-1999_fast_lpjgforcing.nc`
were checksum- and variable-verified as exact copies of their intended CERES
source files.

CERES inputs:

```
/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_200003-202106.nc
/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc
```

Full construction and plotting provenance:

```
/work/ab0995/a270270/analysis/Offline_LPJ_GUESS_analysis/
```

Prepared but not yet run AMIP NH temperature perturbations:

```
/work/ab0995/a270270/input/cruncep_v7/AMIP_noLPJG_1d_1870-1879_TCO95_PI_fluxfix_NHplus1K.nc
/work/ab0995/a270270/input/cruncep_v7/AMIP_noLPJG_1d_1870-1879_TCO95_PI_fluxfix_NHplus2K.nc
```

Only `tas`, `tasmin`, and `tasmax` at latitude >=45 degrees N were incremented;
all other values are unchanged from the source AMIP forcing.

## Observations
```
ERA5 (2m T, seasonal):  /work/ab0246/a270092/obs/era5/netcdf/T2M_{DJF,JJA,DJFM,...}.nc
```

## Static input
```
FESOM CORE3 mesh:       /work/ab0246/a270092/input/fesom2/core3/
LPJG oifs_forcing dir:  /work/ab0246/a270092/input/lpj-guess/oifs_forcing/
OIFS cmip6 forcing:     /work/ab0246/a270092/input/oifs-48r1/cmip6-data/
```

## Model codes & tooling
```
awiesm3-develop:        /work/ab0246/a270092/model_codes/awiesm3-develop
oifsamip-cy48:          /work/ab0246/a270092/model_codes/oifsamip-cy48
lpjg-spinup-develop:    /work/ab0246/a270092/model_codes/lpjg-spinup-develop   (LPJG @ tag 4.1.6 / 2e9ce65)
esm_tools (JanStreffing fork): ~/esm_tools
reval (release_evaluation_tool2): /work/ab0246/a270092/software/release_evaluation_tool2
  key tools: part2_rad_balance.py, part8_t2m_vs_era5.py, part4_cmpi.py
  lpjg_helpers.py (used by the PFT/decisive figure scripts): scripts/lpjg_helpers.py (copied into scripts/figures/)
```

## Original working directory (source of everything here)
```
/work/bb1469/a270092/eval/
  report.tex report.pdf  plots/  scripts/  notebooks/  data/  configs/  run_figures/<run>/
```

## Key parameters referenced in the report
- LPJG establishment gates: `twmin_est=5`, `gdd5min_est=500`, evergreen `tcmin_est=-30` (`global.ins`).
- LPJG competition knobs (Part II verdict → the lever): boreal `greff_min` (0.03 BNE / 0.09 BINE,BNS),
  C3-grass competitiveness, cold-tree `pstemp_low=10`.
- OpenIFS stable-BL (Phase 1B): `ZLMIN` (`vdfexcu.F90:201`), `ZCB` (`vdfexcu.F90:190`).
- XIOS flux de-accumulation period fix: `field_def_lpjg_safe.xml.j2`, divisor `= NFRHIS*TSTEP` (report App. B.2).
- Follow-up coupled branch candidates (2026-07-24): 06T (`oasis3mct.time_step=3600`, `use_momix=True`, `momix_kv=0.01`, `kpp_av0=0.003`, `kpp_kv0=0.003`, `kpp_kvbckg=5.0e-06`) and 06V (06T + `ENTSTPC3=1`), then LPJ levers 07A (`BNS greff_min=0.035`) and 07C (`C3G pstemp_low=12`).
