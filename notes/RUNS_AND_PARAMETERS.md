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

## 1. THE KEPT CONFIGURATION — G4

**Carry this forward.** Best configuration of the campaign on every headline metric,
seasonally clean, and it repairs most of the Nordic Seas damage G1 introduced.
G4 $=$ G1 $+$ tundra, i.e. F4 (`RVRSMIN` 250→1000 for veg types 3, 4) $+$ D2b
$+$ `RVRSMIN(9)` 80→225.

| metric | control | G1 | **G4** | Δ vs control | target |
|---|---:|---:|---:|---:|---:|
| **Siberia JJA T2m** [°C] | 9.73 | 10.25 | **10.68** | **+0.952** (t=7.68) | ≈12.2 |
| Siberia sfc net SW [W/m²] | 153.78 | 159.31 | **161.88** | +8.098 | 166.26 |
| Siberia cloud area [%] | 78.14 | 75.45 | 74.77 | −3.369 | 69.59 |
| **SO SW RMSE** (priority) | 6.877 | 4.809 | **4.800** | −2.076 | — |
| SO TOA SW CRE [W/m²] | −60.29 | −63.13 | −63.04 | −2.753 | −68.14 |
| SO cloud area [%] | 83.07 | 83.59 | 83.55 | +0.484 | 89.72 |
| global net TOA [W/m²] | +0.64 | +0.45 | +0.52 | −0.128 | ~0 |
| tropics net TOA [W/m²] | 42.61 | 42.71 | 42.73 | +0.124 | 45.11 |
| subpolar N Atl SW RMSE | 5.007 | 4.872 | **4.738** | −0.269 | — |
| Nordic Seas SW RMSE | 9.058 | 9.358 ⚠ | **9.147** | **+0.089** | — |
| global T2m RMSE [K] | 1.579 | 1.553 | **1.543** | −0.036 | — |
| DJF / MAM / SON T2m | — | +0.07/−0.30/+0.28 | −0.385 / −0.360 / +0.185 | all **within** noise | |

**G4 closes 48 % of the boreal bias**, up from G1's 26 %.

### How to set it — namelist only, no source edits

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
                "ECE_TUNE_RVRSMIN(9)": 225.0    # tundra, 80 -> 225 (table consistency)
```

⚠ **`NAMSURFTUNE`, not `NAMECECFG`.** `NAMECECFG` is read **twice** from the same
`fort.4` — by `ECE_CONFIG` in `arpifs/ecearth.F90` and by `SURFECE_CONFIG` in `surf` —
and a Fortran namelist read aborts with *invalid reference to variable* on any name the
reading module does not declare. Putting the tune entries in `NAMECECFG` kills the arpifs
read at `su0yoma.F90:152`. Verified the hard way by a 1-day test run.

Runscript: `oifsamip-cy48-levante-TCO95L91_G4_tundra.yaml`.

### The three components

| lever | change | mechanism | boreal JJA | SO SW RMSE |
|---|---|---|---:|---:|
| **F4** | `RVRSMIN` 250→1000, types 3, 4 | less transpiration → more sensible heat | +0.749 (t=6.06) | −0.217 |
| **D2b** | `RCL_INPSEA=0.2`, `RCL_INPPMIN=700 hPa` | fewer ocean ice nuclei → more supercooled liquid → brighter cloud | −0.222 | **−1.914** |
| **tundra** | `RVRSMIN(9)` 80→225 | same as F4, on the 25.6 % of the box F4 misses | **+0.431** (G4−G1) | −0.009 |

**Why tundra is justified — a table-consistency argument, independent of our bias.**
`RVRSMIN(9) = 80 s/m` is the **lowest value of any vegetated type in HTESSEL**: below
crops, short grass and tall grass (all 100), semidesert 150, deciduous broadleaf 175,
evergreen/deciduous shrubs 225, bogs/marshes 240, desert and needleleaf 250. The model
therefore has arctic tundra transpiring more freely than tropical grassland. Its closest
physiognomic analogues in the same table are the shrubs at 225, and bogs/marshes — which
co-occur with tundra in this very box — are 240. Contrast F4, whose 250→1000 remains an
unanchored 4× excursion.

Vegetation cover in the box, area-weighted over land (model's own `tvh/tvl/cvh/cvl`):
**tundra 25.6 %**, decid. needleleaf 19.2 %, mixed forest 5.2 %, bogs/marshes 3.9 %,
evergr. needleleaf 3.4 %, evergr. shrubs 3.0 %. F4 reaches ~24 %; tundra is another 26 %.

### ⚠ Known costs, caveats, and one falsified mechanism

- **The melt mechanism claimed for G4 is FALSIFIED.** The prediction on record was a
  measurable drop in May–June snow water equivalent. Measured (thresholds ±5.25 May,
  ±2.72 Jun): **G4 May +2.75, Jun +1.60 — both within noise.** T2m rose with the snow
  mass unchanged, so G4's gain is the **plain sensible-heat route**, the same as F4.
  The lever works and is independently justified, but *not for the reason given*.
- **G1/F4 significantly *increases* May–June snow** (+5.58, +5.68, both clearing
  threshold). F4 improves T2m while **aggravating** the snow bias — plausibly by cutting
  evapotranspiration and therefore sublimation, a real snowpack sink. The warming partly
  offsets a bias it also worsens. Unresolved.
- **`RVRSMIN` does not saturate — do NOT go past 1000.** G2/G1/G3 at 500/1000/2000 give
  +0.336/+0.521/+0.876: increments per doubling are +0.185 then **+0.355**, i.e.
  *accelerating*, no knee. 1000 is defensible only as *before the winter damage starts*,
  not as "on the knee". G3 at 2000 costs **DJF −0.798**, clearing ±0.588 — genuine winter
  damage, joining B5 and H2.
- SO **cloud-area** deficit (~6 pp) untouched by everything ever tried.
- Tropics net TOA still ~2.4 W/m² low.

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

### Round 14 results (2026-08-03) — G4 adopted, F4 magnitude capped, mechanism falsified

All three were **namelist-only**, sharing one binary with no rebuild — the first round the
`NAMSURFTUNE` work made possible.

| run | setting | Siberia JJA | DJF | verdict |
|---|---|---:|---:|---|
| **G2** | `RVRSMIN(3,4)=500` | +0.336* | −0.325 | half-strength; 64 % of G1's gain |
| **G3** | `RVRSMIN(3,4)=2000` | +0.876* | **−0.798\*** | **rejected** — winter damage |
| **G4** | G1 + `RVRSMIN(9)=225` | **+0.952\*** (t=7.68) | −0.385 | **ADOPTED** |

**Predictions vs outcome, both recorded before the runs finished:**

| prediction | outcome |
|---|---|
| G2 keeps >60 % of G1's gain if saturating | ✅ 64 % |
| G3 adds < +0.15 K beyond G1 if saturating | ❌ **added +0.355** — not saturating |
| G4 gains +0.2…+0.5 K beyond G1 | ✅ +0.431 |
| G4 gain accompanied by a **drop in May–June SWE** | ❌ **falsified** — May +2.75, Jun +1.60, both within noise |

So the boreal gain is real and reproducible, but **the spring-melt mechanism that motivated
G4 is not supported**: the temperature rises without the snow moving. Two mechanisms for the
June albedo bias have now been falsified — snow *cover formulation* (round 13) and snow
*melt rate* (round 14).

### Not tuning levers — never add to `RUNS`

`amip_A2_kkland150` (superseded), and 11 LPJG forcing-generator runs from 3–5 July 2026
(`amip_lpjgforce_chk`, `amip_nolpjg_forc`, `amip_nolpjg_forcing`, `amip_nolpjg_pi1870`,
`amip_pi_clean1/2`, `amip_pi_dbg1/2/3`, `amip_pi_fixtest`, `amip_pi_forcing`). These
emit the daily LPJG forcing set, not evaluation fields, so they legitimately have no
`atm_remapped_1m_2t_*`. Listed in `NOT_LEVERS` in `runs.py`.

**Also: do not touch `fgdbg02`.**

---

## 3b. The causal web, and what is still dark

Status markers: **[OBS]** established against observation (CERES, Rutgers/IMS satellite);
**[MOD]** model-vs-model only (ERA5 is an HTESSEL sibling — suggestive, not authoritative);
**[DEAD]** tested and falsified; **[?]** unknown.

```
        ROUTE A  (transpiration -> cloud)          ROUTE B  (snow -> albedo)
        =================================          =========================
 RVRSMIN too low (tundra 80 s/m)             [?] unknown Mar-Apr energy sink
        |                                             |
        v                                             v
 over-transpiration                    [MOD]   melt-out 13 days late          [OBS]
   +-> latent cooling                                 |
   +-> moist boundary layer                           v
        |                                      Jun snow cover 2x observed      [OBS]
        v                                             |
 JJA cloud +8.5 pp                     [OBS]          v
        |                                      Jun albedo +0.046               [OBS]
        v                                             |
 -7.0 W/m2 surface SW                  [OBS]   -5.6 W/m2 surface SW            [OBS]
        |                                             |  [?] does this reach T2m?
        +------> Siberian JJA T2m ~2.0 K too cold <---+
                          |
                          |  (severed deliberately in AMIP)
                          v
              forest fails -> tundra -> brighter+smoother -> colder  (coupled runaway)

 G4 CLOSES ROUTE A: cloud term -7.04 -> +1.20, T2m +0.952 K, 48 % of the bias
 ROUTE B UNTOUCHED: albedo term -5.45 -> -5.59, and G4 makes the snow WORSE
```

**Route A is closed.** Excess stomatal conductance cools latently and moistens the
boundary layer → low cloud → SW removed → cold. G4 removes the cloud term entirely.
Defensible on table-internal grounds for tundra; **not** for needleleaf, where
`RVRSMIN` 250→1000 is unanchored and does not saturate.

**Route B is open**, and everything easy is already dead: not the cover formulation
(round 13), not energy supply (model absorbs *more* net SW than CERES in Mar–Apr), not
snow brightness (model snow is *darker* than observed), probably not cold content [MOD].
**[?] Where the March–April energy goes is the central unknown.** Candidates: turbulent
flux (unconstrained by any observation), conduction/refreeze in the pack, or the
canopy-shaded snow tile (7).

**⚠ The link that may not exist.** G4 raised T2m by +0.952 K while May–June snow mass
did **not** change and snow cover and melt date got **worse**. Temperature moved a full
kelvin with the snow term stationary or degrading. Either the June albedo error is
weakly coupled to T2m, or the routes offset in a way we have not separated. **Until this
is resolved the prize for fixing Route B is unknown** — it may buy far less than its
−5.6 W/m² suggests.

**⚠ New coupled risk.** G4 delays melt-out 24 → 28 May. Later melt is later
growing-season onset for LPJ-GUESS, so G4 may improve AMIP T2m while *harming* forest
establishment through a pathway AMIP cannot see. Argues for running coupled with G4
sooner rather than later.

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
   realistic would make June *worse*. ~~The lever should be snow mass / melt rate.~~
   **← falsified by round 14, see item 2.**
   (Caveat: ERA5's albedo scheme is a relative of HTESSEL, so `asn` agreement is not
   independent evidence; its snow *extent* is observationally constrained via IMS.)

   **Narrowed further (round 14).** Accumulation is right — through the accumulation
   season the model is 7–9 % high and monthly increments track ERA5 almost exactly.
   The divergence is spring: **April the pack still gains (+9.4 mm) while ERA5 has
   peaked (−1.6)**, and **May melts 30 % too slowly** (−33.5 vs −47.5). Against CERES,
   May's SW↓ is nearly right (−3.5) while SW_net is **−13.9** → May is almost pure
   surface albedo, not cloud. Self-reinforcing: late melt → more May snow → brighter →
   less absorbed → less melt energy → snow into June.

2. **⚠ The spring-melt mechanism is ALSO falsified (round 14).** G4 raised Siberian JJA
   T2m by +0.952 K with **no** significant change in May–June snow water equivalent
   (May +2.75 against ±5.25, Jun +1.60 against ±2.72). Temperature moved; snow did not.
   Two mechanisms are now dead — snow *cover formulation* (round 13) and snow *melt rate*
   (round 14) — and the June albedo bias remains unexplained.

   Worse, **G1/F4 significantly *increases* May–June snow** (+5.58, +5.68, both clearing
   threshold), so the campaign's best boreal lever aggravates the snow bias while improving
   temperature. Leading suspect: cutting evapotranspiration also cuts **sublimation**, a
   real snowpack sink. Untested.

   Remaining candidates, none yet tried: the **snow-albedo decay timescale**; the
   **canopy-masking (tile 7) formulation**, which the round-13 residual showed is already
   doing most of the work in the melt months; and the **`RVVEGALB` vegetation albedo
   table**, which the campaign has never touched.
3. **SO cloud-area deficit ~6 pp** — untouched by every lever tried.
4. **The energy target, now decomposed PER RUN (2026-08-03).** *Where* the imbalance lives
   was settled long ago — the global budget is right **by cancellation**, with the Southern
   Ocean carrying +5.8 W/m² and contributing +0.57 on its own, from a cloud deficit (the
   Hyder signature). That is unchanged. What was missing: the decomposition existed **only
   for the control**, so across 33 runs the energy target was a single global number with no
   idea which band a lever moved. Now in the standard table via `evaluate.sh`.

   Δ net TOA vs control [W/m²]:

   | | 60–90S | 30–60S | tropics | 30–60N | global |
   |---|---:|---:|---:|---:|---:|
   | *gap to close (CERES−ctl)* | **+8.38** | −0.61 | −2.50 | +2.51 | |
   | A1a ovl=0.10 | −3.33 | −2.80 | −1.06 | −2.44 | −1.90 |
   | **D2a** inpsea.2 | **−2.83** | −1.91 | −0.38 | −0.79 | −0.93 |
   | D2b inp+p700 | −1.76 | −0.83 | +0.09 | −0.19 | −0.28 |
   | B3 clddiff | −0.11 | +0.89 | **+2.02** | +0.86 | +1.34 |
   | B5 capdcycl0 | −0.47 | +0.06 | **−1.93** | −0.14 | −1.00 |
   | F4 rsmin1000 | −0.40 | +0.04 | +0.07 | +0.30 | +0.10 |
   | **G4** | −1.98 | −0.65 | +0.12 | +0.28 | −0.13 |

   - **60–90S is the most responsive band in the model**, moving −0.4 to −3.3 under levers
     designed for other purposes. Not inert — just never targeted.
   - **G4 already buys ~−2.0 of the +8.4**, all from D2b; roughly a quarter of the SO TOA
     gap closed as a by-product. The rest is the same ~6 pp cloud-area deficit.
   - **Surface levers are regionally clean**: F4 +0.07 and G4 +0.12 in the tropics, against
     B3's +2.02 and B5's −1.93. The tropical guardrail is now visible as a band.

   ⚠ **The PI target is net TOA ≈ 0, NOT the CERES column.** CERES is present-day and carries
   the real warming imbalance (~+1.1 W/m² globally), so the control's +0.64 is *too positive*
   even though it sits *below* CERES. Use per-band CERES differences to locate error, never
   as a global tuning target.
5. **Tropics net TOA 2.4 W/m² too low** (42.6 vs 45.1).
6. **Nordic Seas SW RMSE** — G1 degraded it by +0.300 against a ±0.052 threshold, the only
   lever ever to move that box. **G4 largely repairs it** (+0.089), unexpectedly and for no
   reason we understand: adding tundra resistance should not act on a Nordic ocean box.
   Worth a look before trusting it — an unexplained improvement is as much a loose end as
   an unexplained degradation.
7. `RVRSMIN` 4× is unanchored **and does not saturate**. G2/G1/G3 at 500/1000/2000 give
   +0.336/+0.521/+0.876 — increments per doubling of +0.185 then **+0.355**, accelerating,
   with no knee to stop at. 1000 is justified only as *the largest value before the winter
   damage appears* (G3 at 2000 costs DJF −0.798, clearing ±0.588). Anchoring this against a
   transpiration observation is now **more** important, not less.
8. Winter surface albedo is ~0.02–0.03 **too low** (opposite sign to June) — a
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
