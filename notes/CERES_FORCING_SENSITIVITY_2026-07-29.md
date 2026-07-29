# CERES-constrained radiation forcing sensitivity — 2026-07-29

## Scientific question

Does the surface-radiation treatment explain the different equilibrium boreal
vegetation obtained with AMIP and CRUNCEP forcing?

## Experiment

Four otherwise identical 2000-year TCO95 offline LPJ-GUESS spin-ups were
compared in model year 3900:

1. OpenIFS–AMIP 1870–1879.
2. Calibrated CRUNCEP v3 1901–1910.
3. CRUNCEP v3 with `rsns` and `rlns` replaced by monthly CERES EBAF Ed4.1
   fields from 2001–2010 and linearly interpolated to daily values.
4. CRUNCEP v3 with monthly radiation means constrained to the official CERES
   climatology while retaining CRUNCEP relative daily shortwave and additive
   daily longwave anomalies.

Temperature, precipitation, humidity, wind, coordinates, and time are unchanged
between CRUNCEP v3 and both CERES-derived forcings.

## CERES provenance

- Product: NASA Langley CERES EBAF Edition 4.1.
- DOI: `10.5067/TERRA-AQUA/CERES/EBAF_L3B004.1`.
- Direct source:
  `/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_200003-202106.nc`.
- Climatology/reference:
  `/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc`.
- Reference variable: `sfc_net_sw_all_clim`, all-sky surface net shortwave.
- Official climatology period: July 2005–June 2015.

The direct-forcing experiment uses monthly CERES values from 2001–2010. The
daily-variability experiment uses the official 12-month climatology.

## Result

Final-year NH-land display-grid mean of the annual `treeFrac` diagnostic:

| Forcing | Tree fraction | Difference from AMIP |
|---|---:|---:|
| AMIP | 0.186 | 0.000 |
| CRUNCEP v3 | 0.269 | +0.083 |
| CRUNCEP + CERES direct | 0.262 | +0.076 |
| CRUNCEP + CERES daily variability | 0.263 | +0.076 |

Replacing CRUNCEP radiation with CERES lowers the mean by only 0.006–0.007 and
leaves about 92% of the CRUNCEP–AMIP contrast. The two temporal treatments are
nearly indistinguishable in the hemispheric mean.

## Interpretation

Surface net radiation is not the dominant source of the AMIP–CRUNCEP
tree-fraction difference. Combined with the controlled full-forcing transfer
test in report Part IV, this locates most of the response in other meteorology
(especially temperature/GDD5), field covariability, and nonlinear LPJG
response. It does not demonstrate a residual competition error. Competition
parameters should be tested only after a forcing-consistent, same-binary
re-spin.

## Important limitations

- CERES is satellite-era, not preindustrial or 1901–1910.
- The local CERES files are monthly; direct daily values are interpolated.
- CERES-derived forcing minus CERES is a construction check, not independent
  validation.
- AMIP, CRUNCEP, and CERES differ in epoch and generation method.

Detailed forcing algorithms, validation ranges, common-grid mapping, masks, and
plots are retained in:

`/work/ab0995/a270270/analysis/Offline_LPJ_GUESS_analysis/`

## Next controlled experiment

Two AMIP forcing files are ready with exactly +1 K and +2 K applied to `tas`,
`tasmin`, and `tasmax` at latitude >=45°N. No LPJ-GUESS spin-up results exist
for these inputs yet. They should be run with all LPJG parameters fixed to test
temperature thresholds separately from radiation and model-parameter changes.
