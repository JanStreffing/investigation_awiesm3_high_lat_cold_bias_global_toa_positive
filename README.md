# AWI-ESM3 high-latitude cold bias & positive global TOA imbalance — investigation

Post-processing, figures, scripts and the write-up for the AWI-ESM3 (v3.4,
CORE3 / TCO95) PI-control investigation.

**The report is `report/report.pdf`** (source `report/report.tex`, recompilable
in place). It is structured in three parts:

1. **Tuning away the high TOA imbalance** — 22 runs; every effective lever cools
   the model, and the aggressive coupling/ocean-mixing levers cut the +1.5 W/m²
   imbalance hardest but push the planet too cold vs ERA5. The aggressive lever is
   traced to a Southern-Ocean open-ocean-convection mechanism.
2. **Counter-tuning the NH cold bias (the structural floor)** — the cold bias that
   makes the aggressive lever expensive is largely *not* the tuning: it is a
   structural LPJ-GUESS boreal-forest collapse. A clean WITH-vs-WITHOUT-LPJG
   experiment attributes it (snow-albedo + roughness→inversion, plus an
   over-strong stable-BL residual). A **decisive feedback-free standalone re-spin**
   then shows the forest deficit is an LPJG **tree-vs-grass competition** problem,
   **not** a climate one (`plots/v3_decisive_climate_vs_veg.png`).
3. **The path forward** — raise the structural warm floor (vegetation competition +
   stable-BL mixing), *then* spend a moderate cooling lever to reach surface
   balance without the RMSE blow-up.

**2026-07-24 follow-up:** a finished 50-year coupled run started from the
CRUNCEP3-initialized vegetation state (alias "080a") refines the picture. The
failure is **not** a simple NH-wide cold drift. By year 50 the
system shows a **warm-ocean / cold-land split**, and the vegetation loss is
strongest over **Siberia / East Siberia**, where `TREEFPC` and `AGDD5` both
decline even though NH45+ aggregate thermal metrics can increase. The current
working interpretation and next-run proposal are recorded in:

- `notes/COUPLED_CRUNCEP3_INIT_BRANCH_FOLLOWUP_2026-07-24.md`

See **`DATA_PATHS.md`** for every experiment and its Levante path.

## Layout

```
report/            report.tex, report.pdf, report_values.tex, run_figures/<run>/
plots/             all figures used in (and beyond) the report (PNG)
scripts/
  analysis/        per-mechanism analysis (DJF drivers, sea-ice trends, LPJG diagnostics, decisive test)
  figures/         plotting notebooks (SO diagnostics, PFT maps, seasonal cycle) + lpjg_helpers.py
  summary/         ensemble summary/scatter, t2m-bias, CMPI breakdown, report-value generation
  sbatch/          SLURM wrappers (LPJG diag, forcing generation, nc conversion)
data/              derived CSVs (master metrics, CMPI by var/region, tree/grass fractions, ...)
configs/           reval per-run evaluation configs
notes/             LPJG_PFT_investigation.md, PROGRESS_REPORT.md, COUPLED_CRUNCEP3_INIT_BRANCH_FOLLOWUP_2026-07-24.md
```

## Which script makes which figure

| Figure (`plots/`) | Script |
|---|---|
| `summary_scatter_official`, `summary_TOA_timeseries`, `summary_scatter_TOA_vs_T2m` | `scripts/summary/plot_summary.py`, `plot_summary2.py` |
| `t2m_bias_*`, `t2m_bias_multipanel` | `scripts/summary/plot_t2m_bias.py` |
| `cmpi_breakdown_by_{variable,region}` | `scripts/summary/cmpi_breakdown.py` |
| per-run `radiation_budget`, `t2m_vs_ERA5`, `cmpi` (in `report/run_figures/`) | `reval` part2 / part8 / part4 (see `configs/`) |
| `mld_/sic_/sit_winter_antarctic`, `temp_wws_section_tuning` | `scripts/figures/run_mld_sic_sit_tuning.py`, `run_temp_wws_*.py` |
| `pft_boreal_compare`, `pft_lai_boreal_diff` | `scripts/figures/run_pft_boreal_compare.py`, `run_pft_lai_diff_maps.py` |
| `baseline_boreal_evolution` | `scripts/figures/run_baseline_evolution.py` |
| `nolpjg_seasonal_cycle`, `nolpjg_maps` | `scripts/figures/run_nolpjg_compare.py` |
| `djf_roughness_inversion_maps` | `scripts/analysis/plot_djf_roughness_maps.py` (+ `djf_roughness_inversion.py`) |
| `djf_roughness_scatter` | `scripts/analysis/correlate_djf_roughness.py` |
| `djf_bias_vs_era5` | `scripts/analysis/djf_bias_vs_era5.py` |
| `djf_second_driver` | `scripts/analysis/djf_second_driver.py` |
| `establishment_climate_check` | `scripts/analysis/establishment_climate_check.py` |
| `seaice_t2m_trend`, `seaice_volume_trend` | `scripts/analysis/seaice_t2m_trend.py`, `seaice_volume_trend.py` |
| **`v3_decisive_climate_vs_veg`** (the decisive test) | `scripts/analysis/v3_decisive.py` |

2026-07-24 coupled follow-up figures copied into this repo under `plots/`
(`treefpc_coupled_cruncep3_init_y{10,50}_minus_y0_nh`,
`t2m_seasonal_bias_coupled_cruncep3_init*_vs_cruncep3_nh`,
`t2m_bias_vs_base06_all_tuned_coupled_nh_10yrseasonal_{DJF,MAM,JJA,SON}`,
`net_radiation_panel_baseline_A_to_V_07A_to_C_HRgoal_LRpi_last30`) are
generated in the working analysis directory:
`/work/ab0995/a270270/analysis/LR_offline_LPJ_GUESS_tunning_new/scripts/`.

## Recompiling the report
```
module load texlive/live2025-gcc-13.3.0
cd report && latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex
```
(`plots/` is symlinked into `report/`; per-run figures are in `report/run_figures/`.)

## Notes
- Scripts hard-code Levante absolute paths (see `DATA_PATHS.md`); they are
  provenance/reproduction records, not a turnkey pipeline.
- `scripts/figures/lpjg_helpers.py` is copied from the `reval` tool; the PFT and
  decisive-test scripts import it.
