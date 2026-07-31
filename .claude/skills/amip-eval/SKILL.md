---
name: amip-eval
description: Standard evaluation protocol for the AWI-ESM3 AMIP tuning campaign - every metric, region, threshold, reference dataset and known trap. Use when evaluating any oifsamip run, comparing tuning levers, deciding whether a result is significant, or reporting campaign results.
---

# AMIP tuning campaign — standard evaluation

Run **every** section below when evaluating a round. Reporting a subset silently drops
guardrails and has repeatedly produced wrong conclusions in this campaign.

## 0. Before believing anything: the two traps that have burned us

**Trap 1 — the noise floor.** At the original 4-year window the 95 % detection threshold on
Siberian JJA T2m was **±0.89 K** and only 1 of 18 runs cleared it. A whole round was spent
promoting a lever (B8, +0.502 K) that was noise, then building two more runs on top of it.
At 44 years the threshold is **±0.240 K**. *Always run `noise_floor.py` and quote `t`, not
just the delta.* "Inert" and "unresolved" are different claims — B5 looked inert at +0.001 K
over 4 years and is the best boreal lever at +0.407 K over 44.

**Trap 2 — reference period.** The model is 1870s; ERA5's local file is **1990–2014** and
CERES EBAF here is the **07/2005–06/2015** climatology. An indirect estimate of the offset
(HadCRUT5 chain) gave +1.12 K and was **wrong by 2.7×**. The *measured* offset from
`amip_presentday` is **+0.42 K** on Siberian JJA T2m and **−0.07 W/m²** on SO SW CRE. Use the
measurement. Prefer `amip_presentday` (1990–2014) vs ERA5/CERES in their own period whenever
the question allows.

## 1. Scripts, in the order to run them

| script | what it gives |
|---|---|
| `scripts/analysis/eval_round10_A.py` | the main table: SO, Siberia, global guardrails, deep-water + global RMSE |
| `scripts/analysis/noise_floor.py` | run × year ANOVA → detection threshold and `t` for Siberian JJA T2m |
| `scripts/analysis/rmse_significance.py` | same ANOVA on per-year spatial SW RMSE (deep-water boxes) |
| `scripts/analysis/vertical_profiles_prep.sh` then `vertical_profiles.py` | model − ERA5 profiles of T, q, RH |

Add new runs to the `RUNS` list in **all three** of the first scripts — they each carry their
own copy. Runs with missing output are skipped with a printed warning, so adding a
still-running experiment is harmless.

## 2. Metrics, regions, and current targets

Evaluated over **1872–1915** (44 yr) against control `amip_pi_base`, unless the run is
`amip_presentday` (1990–2014).

| metric | region | control | target | note |
|---|---|---:|---:|---|
| TOA SW CRE | SO 45–65S **ocean** | −60.29 | −68.14 (CERES) | gap −7.85 |
| cloud area | SO 45–65S ocean | 83.07 % | 89.72 % | ~6 pp deficit, untouched by everything tried |
| T2m | Siberia 55–75N, 60–180E **land**, JJA | 9.73 °C | ≈12.2 (ERA5) | bias ≈ −2.0 K period-clean |
| surface net SW | same Siberia box | 153.78 | 166.26 (CERES) | |
| **global net TOA** | global | +0.64 | **~0**, but piControl baseline is **+0.79** | |
| global surface flux | global | +0.49 | 0 | includes snow enthalpy `sf × 3.3355e8` |
| tropics net TOA | 30S–30N | 42.61 | 45.11 (CERES) | **model is too LOW — 2.5 deficit** |
| SW RMSE | SO 45–65S ocean | 6.88 | — | **the priority metric** |
| SW RMSE | subpolar N Atl 50–65N, 60W–0, ocean | 5.01 | — | |
| SW RMSE | Nordic Seas 65–80N, 20W–20E, ocean | 9.06 | — | nothing significant, ever |
| global SW RMSE / T2m RMSE | global | 14.20 / 1.58 | — | rank only, not skill |

"Deep-water formation SW RMSE" is a **surface shortwave error over ocean points** in three
boxes — nothing oceanic is simulated in AMIP. The name describes why the boxes were chosen.

Detection thresholds at 44 yr: Siberian JJA T2m **±0.240 K**; SO SW RMSE **±0.019**; subpolar
N Atl **±0.024**; Nordic **±0.053**; global SW **±0.023**.

## 3. Known-broken and known-tricky data

- **Model pressure-level `r` (RH) is identically zero in every run.** Compute RH from `t` and
  `q` instead (Bolton), applied identically to both datasets.
- **No cloud on pressure levels** (`cc`/`clwc`/`ciwc` absent) — RH is the only vertical cloud proxy.
- **`NODE.001_01` is ~3 GB and detected as binary** — plain `grep` silently returns nothing.
  Use `LC_ALL=C grep -a`.
- Model **monthly** `pl` output exists (36 MB/yr) alongside 6-hourly (4.6 GB). Use monthly.
- ERA5 monthly pressure levels are already on DKRZ: `/pool/data/ERA5/E5/pl/an/1M/<param>/`
  (130 T, 133 q, 157 RH); surface: `sf/an/1M/` (167 T2m), fluxes `sf/fc/1M/` (178/179/208/209).
  **No CDS download needed.** ERA5 is reduced-Gaussian: `-setgridtype,regular` first, and
  `-f nc` or cdo writes GRIB. All 19 model levels exist exactly in ERA5's 37.
- ERA5 accumulated fluxes are J/m² — divide by 86400 for W/m².

## 4. Reporting rules that have earned their place

- Quote **`t` or "within noise"** beside every boreal number.
- Report **absolute values**, not only deltas — the reduction 6.88 → 4.96 reads very
  differently from "−0.17".
- Check **tropics net TOA** on any lever touching convection or cloud. B5 wins the boreal at
  a cost of −1.94 W/m² there.
- **Never predict superposition.** AB and ABB8 both got it wrong *in sign*. Measure it.
- State a **falsifiable prediction** before a run and report the outcome against it. D1's
  prediction failed cleanly and that was worth more than a vague success.
- When a claim is retracted, **strike it in place with the reason** rather than deleting.

## 5. Build discipline for source-edit runs

One model tree, so source experiments **serialise**: edit → `esm_master recomp-oifsamip-cy48/oifs-48r1`
→ **verify the library md5 changed** → submit → confirm the experiment staged that md5 → revert
→ next. `comp-` once silently reused a stale object and cost five wrongly-killed runs.

- Detach long builds (`setsid nohup`) so a session teardown cannot leave `install/` empty,
  and add a separate tracked watcher if a notification is wanted.
- Existing experiments keep their own `bin/` and `lib/`, so continuation legs are immune to
  a rebuild; **new** experiments stage from the shared tree.
- Leave the tree at as-released defaults; make new namelist parameters **no-ops by default**.
- `oifsamip` runs use `model_dir = model_codes/oifsamip-cy48`. `recomp-awiesm3-develop-is` is
  the *coupled* tree and does not affect AMIP runs.
- Legs must be ≤ 8 h wallclock. Throughput varies 3–4× with system load; size legs against
  the **slowest** observed rate (~12,700 steps/h → 8-year legs), not a fast night-time sample.

## 6. Forcing configuration — verify, don't assume

`NCMIPFIXYR` must appear in **`&NAMECECMIP`**. 48r1 renamed the old `&NAMECECMIP6`, which
survives in the esm-tools fort.4 template, so patching it is silently ignored and the run
falls back to `NCMIPFIXYR = -1`, i.e. **transient** forcing. Set `oifs.scenario: "piControl"`
and check the generated `fort.4`, plus `CO2 =` in `NODE.001_01` (should be constant ≈284.3).
