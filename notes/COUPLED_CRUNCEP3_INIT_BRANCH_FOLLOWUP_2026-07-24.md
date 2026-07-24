# Coupled CRUNCEP3-initialized branch follow-up (2026-07-24)

## Purpose

This note records the current understanding after the coupled
CRUNCEP3-initialized branch (alias "080a") completed **50 years**.

The goal of this follow-up was to test whether the cold-bias / boreal-forest
problem remains simply an inherited LPJ-GUESS spin-up issue, or whether the
fully coupled adjustment adds a new regional climate failure.

## Main result

The finished 50-year coupled run shows that the problem is **not** a simple
NH-wide cooling and is **not** adequately captured by a pure
"competition-not-climate" story.

The clearest pattern is:

- **warmer ocean**,
- **colder boreal land**,
- and persistent **tree-fraction collapse over Siberia / East Siberia**.

So the current best diagnosis is:

1. the LPJG spin-up / competition structure still matters,
2. but the coupled model also develops a **regional land-cooling response** that
   amplifies the vegetation loss.

## Year-50 evidence used

Diagnostics were built from:

- LPJ shard outputs for years `1350`, `1360`, `1399`
- remapped OIFS monthly `2t` fields for `1350-1359` and `1390-1399`
- CRUNCEP3 forcing climatology

Key diagnostic figures (now copied into this repo under `plots/`):

- `treefpc_coupled_cruncep3_init_y50_minus_y0_nh.png`
- `t2m_seasonal_bias_coupled_cruncep3_init_1390_1399_vs_cruncep3_nh.png`

Provenance source directory:

- `/work/ab0995/a270270/analysis/LR_offline_LPJ_GUESS_tunning_new/plots/diff_from_amip_and_era5/`

## Quantitative summary

### NH45+

- `TREEFPC`: `0.336 -> 0.215` (`d = -0.121`)
- `GRASSFPC`: `0.174 -> 0.226`
- `AGDD5`: `1049.9 -> 1132.3`

Interpretation:

- NH45+ aggregates alone are misleading.
- The hemisphere as a whole is not simply becoming colder in all relevant
  thermal metrics.

### Siberia (`55-75N, 60-180E`)

- `TREEFPC`: `0.310 -> 0.121` (`d = -0.188`)
- `Natural_sum`: `0.630 -> 0.425`
- `Total`: `0.604 -> 0.408`
- `AGDD5`: `739.8 -> 650.9` (`d = -89.0`)

### East Siberia (`55-75N, 90-160E`)

- `TREEFPC`: `0.307 -> 0.075` (`d = -0.232`)
- `Natural_sum`: `0.626 -> 0.363`
- `Total`: `0.616 -> 0.357`
- `AGDD5`: `724.8 -> 591.4` (`d = -133.4`)

Interpretation:

- The strongest failure is **East Siberia**, not the entire NH.
- The vegetation loss is accompanied by a real loss of local thermal growing
  potential there.

## What the seasonal temperature plot adds

The final-decade (`1390-1399`) seasonal T2m bias map versus CRUNCEP3 forcing
shows a **warm-ocean / cold-land** split.

Important caveat:

- the ocean portion of that plot should not be used as a strict objective,
  because it compares a coupled atmosphere field over ocean and land against a
  CRUNCEP3 forcing reference field.
- But the **land** signal is still useful and consistent with the boreal
  vegetation decline.

Working interpretation:

- this is not just an inherited spin-up state,
- the coupled system develops a regional continental-land response that remains
  too cold over Siberia by year 50.

## Updated tuning implication

The next tuning path should **not** start with LPJ-only retuning.

Instead, use a two-stage path:

### Stage 1: choose the better physical coupled branch

Test:

- `06T`
- `06V`

Reason:

- these are among the best existing branches for reducing excess global net
  radiation while staying near the desired target.
- `06V` is the preferred candidate; `06T` is the cleaner fallback if the
  `ENTSTPC3=1` atmosphere change proves regionally harmful.

### Stage 2: once the physical branch is chosen, add vegetation support

Then test:

- `07A` first (`BNS greff_min = 0.035`)
- then `07A + 07C` (`C3G pstemp_low = 12`)

Reason:

- `07A` is the most direct boreal-tree persistence lever.
- `07C` is the next-best support lever because it reduces cool-grass pressure.

## Recommended next-run matrix

1. `06T-only`
2. `06V-only`
3. `best(06T/06V) + 07A`
4. `best(06T/06V) + 07A + 07C`

## Acceptance criteria

Only keep a branch if it improves both:

1. **global radiation balance**
2. **boreal land climate / forest persistence**

Specifically check:

- year-50 `TREEFPC` decline in Siberia / East Siberia
- East Siberian `AGDD5`
- MAM and JJA land temperature bias over boreal land
- global net-radiation last-30-year mean

## Bottom line

The 50-year coupled follow-up means the investigation should now be stated more
carefully:

- the old LPJG structural explanation is still part of the story,
- but by year 50 there is also a **regional coupled land-climate problem**,
- so the credible path forward is **physical branch first, LPJ compensation
  second**.