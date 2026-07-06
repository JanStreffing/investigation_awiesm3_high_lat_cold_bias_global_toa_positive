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
