# AWI-ESM3 v3.4 AMIP tuning — runs and kept parameters

Companion to `ATMOSPHERE_TUNING_LOGBOOK.md` (which is chronological). This file is
the **reference**: what every run was, what it measured, and — most importantly —
**which parameter settings we keep** and how to set them.

Setup: `oifsamip` = OpenIFS 48r1 + XIOS + OASIS, TCO95L91, no FESOM. Evaluated over
**1872–1915** (44 yr) against control `amip_pi_base`.
All runs are 46–50 yr; `amip_presentday` is 1990–2014 and is the period-clean pair
for ERA5/CERES comparisons.

Detection thresholds at 44 yr — **quote `t`, never a bare delta**:

| metric | 95 % threshold |
|---|---|
| Siberian **JJA** T2m | ±0.242 K |
| Siberian **DJF** T2m | ±0.588 K |
| Siberian **MAM** T2m | ±0.386 K |
| Siberian **SON** T2m | ±0.431 K |
| SO SW RMSE | ±0.019 |
| subpolar N Atl SW RMSE | ±0.024 |
| Nordic Seas SW RMSE | ±0.052 |
| global SW RMSE | ±0.023 |

---

## 1. THE KEPT CONFIGURATION — G1

**Carry this forward.** Best configuration of the campaign; the only combination that
improves both targets, and the only one that superposes additively.

| metric | control | **G1** | Δ | target |
|---|---:|---:|---:|---:|
| **SO SW RMSE** (priority) | 6.877 | **4.809** | −2.067 | — |
| SO TOA SW CRE [W/m²] | −60.29 | −63.13 | −2.84 | −68.14 |
| **Siberia JJA T2m** [°C] | 9.73 | **10.25** | **+0.521** (t=4.22) | ≈12.2 |
| Siberia sfc net SW [W/m²] | 153.78 | 159.31 | +5.54 | 166.26 |
| global net TOA [W/m²] | +0.64 | +0.45 | −0.195 | ~0 |
| tropics net TOA [W/m²] | 42.61 | 42.71 | +0.096 | 45.11 |
| global T2m RMSE [K] | 1.579 | 1.553 | −0.026 | — |
| Nordic Seas SW RMSE | 9.058 | 9.358 | **+0.300** ⚠ | — |
| DJF / MAM / SON T2m | — | — | +0.070 / −0.296 / +0.280 | all within noise |

### How to set it — namelist only, no source edits (since 2026-08-02)

```yaml
oifs:
    add_namelist_changes:
        fort.4:
            NAMCLDP:
                RVICE: 0.16              # match the coupled tuning runs explicitly
                RCL_INPSEA: 0.2          # D2b: ocean ice nuclei -> 20 %
                RCL_INPPMIN: 70000.0     # D2b: only below ~700 hPa
            NAMSURFTUNE:
                "ECE_TUNE_RVRSMIN(3)": 1000.0   # F4: evergreen needleleaf
                "ECE_TUNE_RVRSMIN(4)": 1000.0   # F4: deciduous needleleaf
```

⚠ **`NAMSURFTUNE`, not `NAMECECFG`.** `NAMECECFG` is read **twice** from the same
`fort.4` — by `ECE_CONFIG` in `arpifs/ecearth.F90` and by `SURFECE_CONFIG` in `surf` —
and a Fortran namelist read aborts with *invalid reference to variable* on any name the
reading module does not declare. Putting the tune entries in `NAMECECFG` kills the arpifs
read at `su0yoma.F90:152`. Verified the hard way by a 1-day test run.

Runscript: `oifsamip-cy48-levante-TCO95L91_G1nml.yaml`.

### The two components

| lever | change | mechanism | boreal | SO SW RMSE |
|---|---|---|---:|---:|
| **F4** | `RVRSMIN` 250→1000 for veg types 3, 4 | less transpiration → more sensible heat | **+0.749** (t=6.06) | −0.217 |
| **D2b** | `RCL_INPSEA=0.2`, `RCL_INPPMIN=700 hPa` | fewer ice nuclei over ocean → more supercooled liquid → brighter cloud | −0.222 | **−1.914** |

**Superposition, measured not assumed** (AB and ABB8 both got it wrong *in sign*):

| | F4 | D2b | predicted | **G1 actual** |
|---|---:|---:|---:|---:|
| Siberia JJA T2m | +0.749 | −0.222 | +0.527 | **+0.521** |
| SO SW RMSE | −0.217 | −1.914 | −2.131 | −2.067 |
| SO TOA SW CRE | −0.079 | −2.640 | −2.719 | −2.840 |
| global net TOA | +0.095 | −0.283 | −0.188 | −0.195 |

Additive to within noise throughout. The plausible reason it holds here and failed
before: F4 is a land surface flux in boreal summer, D2b is ice nucleation over ocean
below 700 hPa. AB/ABB8 combined levers touching the *same* cloud scheme in the *same*
regime.

### Known costs and caveats of G1

- **Nordic Seas SW RMSE +0.300** against a ±0.052 threshold, and *worse* than additive
  (+0.184 expected). No other lever in 31 runs has moved this box at all. Small box,
  not the deep-water priority, but do not let it grow.
- `RVRSMIN` 250→1000 is an **unanchored 4× excursion** from what ECMWF ships. It is
  defensible as a mechanism (see the F4 cause-vs-symptom tests) but it is not
  anchored in a measurement.
- G1 closes **26 %** of the boreal bias and **36 %** of the SO CRE gap. The SO
  **cloud-area** deficit (~6 pp) is untouched by everything ever tried.

---

## 2. Namelist-settable tuning (added 2026-08-02)

The ECMWF `surf` library deliberately ships no namelist, so every HTESSEL lever used
to be a source edit plus a rebuild — which serialises experiments against the single
model tree and puts AWI tuning inside upstream files where it collides with EC-Earth.

`surf/module/surfece.F90` is an AWI module *inside* surf that already carried
`ECE_CPL_LPJG`, the LANDICE API and its own `NAMECECFG` read. The tuning-table overrides
live there too but in their **own `&NAMSURFTUNE` group** (see the warning above for why
they cannot share `NAMECECFG`), applied by `SURFECE_APPLY_TUNING` from `susurf.F90`
**after** `SUSURF_CTL` has filled the tables. The `&NAMSURFTUNE` read is *optional* — a
`fort.4` without the block is not an error, the sentinel defaults simply leave the tables
as released.

| namelist variable | overrides | index | former source lever |
|---|---|---|---|
| `ECE_TUNE_RVRSMIN(0:20)` | min stomatal resistance [s/m] | veg type | **F4** |
| `ECE_TUNE_RVLAI(0:20)` | leaf area index | veg type | F2 |
| `ECE_TUNE_RVCOV(0:20)` | vegetation cover fraction | veg type | F3 |
| `ECE_TUNE_RVZ0H(0:20)` | roughness length for heat [m] | veg type | F1 |
| `ECE_TUNE_RVLAMSK(0:20)` | unstable skin conductivity | veg type | B8 / E1 |
| `ECE_TUNE_RQSNCR` | inverse critical snow depth [1/cm] | scalar | H (rejected) |

Design points that matter:

- Every element defaults to a sentinel (`-999`) meaning **keep the as-released value**.
  A run setting none of them is bit-identical to the untouched code.
- Only **one line** was added to an ECMWF file (`susurf.F90`), next to the existing
  AWI `CALL SURFECE_CONFIG`. Everything else is inside AWI-owned `surfece.F90`.
  This is what makes the tuning mergeable with EC-Earth instead of conflicting.
- Every applied override is written to `NULOUT`, because a silent parameter change is
  exactly what later gets attributed to the wrong lever.
- Nothing in the setup path derives from these entries, so a late override is
  equivalent to editing the table in place. Verified by reading every use of
  `RVRSMIN` in the tree.

HTESSEL vegetation type indices: 1 crops, 2 short grass, **3 evergreen needleleaf**,
**4 deciduous needleleaf**, 5 deciduous broadleaf, 6 evergreen broadleaf, 7 tall grass,
8 desert, 9 tundra, 10 irrigated crops, 11 semidesert, 12 ice caps, 13 bogs/marshes,
14 inland water, 15 ocean, 16 evergreen shrubs, 17 deciduous shrubs, 18 mixed forest,
19 interrupted forest, 20 water/land mixture.

### Still requiring a source build

| change | file | status |
|---|---|---|
| `RCL_INPSEA` / `RCL_INPPMIN` gating | `cloudsc.F90`, `yoecldp.F90`, `sucldp.F90`, `namcldp.nam.h` | **committed**, no-op at defaults, namelist-settable |
| `RCAPDCYCL` mode 4 | `cumastrn.F90` | **committed**, no-op at defaults (D1 falsified) |

---

## 3. Every run in the campaign

`scripts/analysis/runs.py` is the single source of truth for the evaluators — add a
run there and every script picks it up.

### Rejected but instructive

| run | change | why rejected |
|---|---|---|
| **B5** `capdcycl0` | `RCAPDCYCL=0` | JJA +0.407 but **DJF −0.720** (clears ±0.588) — the only lever with genuine winter damage; also tropics −1.94 |
| **B8** `lamsk5` | `RVLAMSK` → 5 | +0.502 at 4 yr looked like the best boreal lever; **at 44 yr it is −0.038, pure noise.** Two runs were built on top of it before the noise floor existed |
| **A1a** `ovl=0.10` | `RCL_OVERLAPLIQICE` 0.10 | buys 82 % of the SO gap at **−0.749 K** in Siberia; a global cloud change, not a regional lever |
| **D1** `capdcycl4` | new land closure | prediction failed cleanly: −0.170, no improvement |
| **H1/H2** `snowcr30` | `RQSNCR` 1/10→1/30 | **round 13, falsified** — see below |
| B3, B4, B6, B7, C1, C2, E1, A1b, A1c, A2, expA, F1, F2, F3 | various | no significant boreal gain, or gain with an unacceptable guardrail cost |
| **F5** `all four` | F1+F2+F3+F4 | +0.746, i.e. **no better than F4 alone** (+0.749) — near-total saturation; F4 carries the whole effect |
| **AB**, **ABB8** | combinations | superposition wrong **in sign** — this is why G1's additivity was measured, not assumed |

### Round 13 (H-series) — the most informative failure

Predicted +0.2 to +0.7 K on Siberian JJA T2m. **Measured +0.020 K, t = 0.17.**
Also **DJF −1.233 K**, the worst winter damage in the campaign.

Why: `ZCVS = min(1, d_cm · RQSNCR)` only responds where snow is **shallow**, and snow
is shallow while *accumulating in autumn*, not while melting in June. The albedo
response landed in Sep–Nov (October −8.2 points) and was −0.95 in June.
**Shallow snow is an autumn phenomenon, not a melt phenomenon.**

The winter cooling is not radiative — Feb/Mar are dark yet cool by −1.65/−1.19 K.
Cutting snow-covered fraction moves area from the snow tiles (5, 7) onto vegetation
tiles, and the snow tile buffers the surface against radiative cooling in polar night.

### Round 14 — in flight (2026-08-03)

All three are **namelist-only**, sharing one binary with no rebuild — the first round the
`NAMSURFTUNE` work made possible.

| run | setting | question | prediction |
|---|---|---|---|
| **G2** | `RVRSMIN(3,4)=500` | is F4 saturating? | if so, keeps >60 % of G1's gain (≥ +0.31 K) |
| **G3** | `RVRSMIN(3,4)=2000` | is F4 saturating? | if so, adds < +0.15 K beyond G1 |
| **G4** | G1 + `RVRSMIN(9)=225` | does tundra close the gap? | +0.2 to +0.5 K beyond G1, **concentrated in June**, with a drop in May/June SWE |

**Falsifiers on record.** G4: if T2m rises with **no** change in May–June snow mass, the
melt mechanism is unsupported and the gain is the plain sensible-heat route. G2/G3: if the
response is **linear** in `RVRSMIN`, we are riding a ramp with no natural stopping point —
an argument *against* the approach, not for a bigger number.

### Not tuning levers — never add to `RUNS`

`amip_A2_kkland150` (superseded), and 11 LPJG forcing-generator runs from 3–5 July 2026
(`amip_lpjgforce_chk`, `amip_nolpjg_forc`, `amip_nolpjg_forcing`, `amip_nolpjg_pi1870`,
`amip_pi_clean1/2`, `amip_pi_dbg1/2/3`, `amip_pi_fixtest`, `amip_pi_forcing`). These
emit the daily LPJG forcing set, not evaluation fields, so they legitimately have no
`atm_remapped_1m_2t_*`. Listed in `NOT_LEVERS` in `runs.py`.

**Also: do not touch `fgdbg02`.**

---

## 4. Open problems

1. **The June surface-albedo bias is still unexplained**, but round 13 plus the ERA5
   decomposition (`scripts/analysis/albedo_decompose.py`) narrowed it sharply. Against
   ERA5 1990–2014, land-masked, Siberian box, June:

   | | model | ERA5 | Δ |
   |---|---:|---:|---:|
   | surface albedo `fal` | 0.2147 | 0.1730 | **+0.042** |
   | snow cover fraction | 0.380 | 0.261 | **+0.119** |
   | snow albedo `asn` | 0.658 | 0.759 | **−0.102** |

   Both f_snow columns use the *same* HTESSEL formula, so the difference is **snow
   amount, not the cover formula**. The model carries roughly twice ERA5's June snow
   mass. And its snow is **too dark**, not too bright — making snow albedo more
   realistic would make June *worse*. **The lever should be snow mass / melt rate.**
   (Caveat: ERA5's albedo scheme is a relative of HTESSEL, so `asn` agreement is not
   independent evidence; its snow *extent* is observationally constrained via IMS.)

   **Narrowed further (round 14).** Accumulation is right — through the accumulation
   season the model is 7–9 % high and monthly increments track ERA5 almost exactly.
   The divergence is spring: **April the pack still gains (+9.4 mm) while ERA5 has
   peaked (−1.6)**, and **May melts 30 % too slowly** (−33.5 vs −47.5). Against CERES,
   May's SW↓ is nearly right (−3.5) while SW_net is **−13.9** → May is almost pure
   surface albedo, not cloud. Self-reinforcing: late melt → more May snow → brighter →
   less absorbed → less melt energy → snow into June.

2. **⭐ Tundra is the largest cover type in the box and no lever has ever touched it.**
   Area-weighted land cover from the model's own `tvh/tvl/cvh/cvl`: **tundra (type 9)
   25.6 %**, deciduous needleleaf 19.2 %, mixed forest 5.2 %, bogs/marshes 3.9 %,
   evergreen needleleaf 3.4 %. F4 reaches ~24 % of the box; tundra is another 26 %.
   `RVRSMIN(9) = 80 s/m` is **the lowest of any vegetated type in HTESSEL** — below
   crops, short grass and tall grass (all 100) — so the model has arctic tundra
   transpiring more freely than tropical grassland, while its shrub analogues are 225
   and co-occurring bogs/marshes are 240. A table-consistency argument, independent of
   our bias. Being tested as **G4**.

3. **SO cloud-area deficit ~6 pp** — untouched by every lever tried.
4. **Tropics net TOA 2.4 W/m² too low** (42.6 vs 45.1).
5. **Nordic Seas SW RMSE** degraded by G1, mechanism unknown.
6. `RVRSMIN` 4× is unanchored; worth trying to anchor against a transpiration
   observation rather than leaving it as a fitted value.
7. Winter surface albedo is ~0.02–0.03 **too low** (opposite sign to June) — a
   separate bias, likely irrelevant to T2m in polar night but real.

---

## 5. Process rules that were learned expensively

- **Noise floor first.** A whole round was spent promoting B8 (+0.502 at 4 yr) which
  is noise at 44 yr. Two further runs were built on it.
- **Use the season's own threshold.** Winter noise is 2.4× summer's. Judging DJF
  against the JJA threshold briefly promoted D2b and B3 to "winter-damaged"; both are
  within the DJF floor.
- **Reference period.** Model is 1870s; ERA5 local file 1990–2014; CERES 07/2005–06/2015.
  The measured offset is **+0.42 K** on Siberian JJA T2m, not the +1.12 an indirect
  HadCRUT5 chain gave (wrong by 2.7×).
- **Never predict superposition — measure it.**
- **State a falsifiable prediction before the run** and report the outcome against it.
  D1 and H1 both failed cleanly and both were worth more than a vague success.
- **Strike retracted claims in place with the reason**, never delete.
- Build discipline: edit → `esm_master recomp-oifsamip-cy48/oifs-48r1` → **verify the
  library md5 changed** → submit → confirm the experiment staged that md5 → revert.
  (`comp-` once silently reused a stale object and cost five wrongly-killed runs.)
  Existing experiments keep their own `bin/`/`lib/`, so continuation legs are immune
  to a rebuild; new experiments stage from the shared tree.
- Legs ≤ 8 h wallclock; throughput varies 3–4× with system load, so size against the
  **slowest** observed rate (~12,700 steps/h → 8-year legs).
