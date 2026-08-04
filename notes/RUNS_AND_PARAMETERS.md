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
| `ECE_SNOW_SCF` | **0** = as-released, **1** = Niu & Yang (2007) depletion, **2** = + SDOR scale-awareness | scalar | round 15 |
| `ECE_SNOW_SCF_Z0` / `_M` / `_RHONEW` | depletion tuning (**0.016 m** / 1.6 / 100 kg m⁻³) | scalar | round 15 |
| `ECE_SNOW_SCF_CSD` | SDOR coefficient, mode 2 only (6.0e-5 m/m) | scalar | round 15 |
| `ECE_LAMSK_SN` | **exposed-snow** skin conductivity `ZSNOW` (7) | scalar | round 16 |
| `ECE_LAMSK_SNGL` / `_SNHV` / `_SNMLT` | `ZSNOW_GLACIER` (8), `ZSNOWHVEG` (20), `ZLARGESN` (50) | scalar | round 16 |

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

### ⭐⭐ Round 15 results (2026-08-04) — the snow fix works, buys no summer, and unmasks a hidden bias

The scheme did **exactly** what it was designed to do:

| Siberian snow cover | satellite | control | G4 | **I1** | **I2** |
|---|---:|---:|---:|---:|---:|
| Apr | 0.950 | 0.962 | 0.961 | 0.931 | 0.931 |
| May | 0.647 | 0.891 | 0.893 | **0.703** | **0.721** |
| Jun | 0.203 | 0.413 | 0.447 | **0.160** | **0.175** |
| **melt-out** | **11 May** | 26 May | 28 May | **12 May** | **13 May** |

Melt-out went from **15 days late to 1 day late**.

**And it delivers no summer temperature.** I2 (snow change alone, no vegetation levers):

| | surface SW gain | JJA T2m | K per W/m² |
|---|---:|---:|---:|
| G4 (vegetation route) | +8.10 | **+0.952** | 0.12 |
| **I2 (snow route alone)** | +2.91 | **+0.094** (noise) | **0.03** |

**`sub:weaklink` is confirmed.** Route B was a genuine, satellite-measured radiative bias
with almost no T2m leverage — 4× less efficient per W/m² than the vegetation route. The
physical reading is clean: energy added over a *melting* snowpack goes into melt, not
sensible heat, which is precisely what it did.

**⚠ The winter response is a COMPENSATING ERROR being removed, not a new one created.**
Scored against ERA5 rather than against the control (period-clean, `amip_presentday`):

| bias vs ERA5 | DJF | MAM | JJA | SON |
|---|---:|---:|---:|---:|
| control | −1.95 | −1.97 | −2.58 | −2.69 |
| G4 | −2.33 | −2.33 | **−1.63** | −2.51 |
| I1 | **−4.71** | −2.28 | **−1.33** | −3.46 |
| I2 | −4.63 | −2.36 | −2.49 | −3.77 |

The control is too cold in **every** season, so there was never a warm bias to protect.

⚠ **RETRACTED (2026-08-04):** the reading that *"excess snow cover was propping DJF up ~2.7 K,
masking a −4.7 K bias"* is **not supported**. Winter cover is *unchanged* between the two
schemes (0.963 vs 0.964, both saturate at 43 cm), so excess **winter** cover cannot have done
the propping. What is measured is only that **I1 cools DJF by −2.76 K by an unidentified
route**. Whether that unmasks a pre-existing bias or creates a new one is **untested** — and
H1 produced the same signature (−1.233 K by cutting snow-tile fraction) and was *rejected*
for it. Identical evidence cannot be disqualifying in one run and revelatory in another.

**Not a radiation-supply problem.** Downward LW over the box, model vs CERES:
Oct **−12.9**, Nov **−11.9**, Dec −3.1, Jan +0.1, Feb +5.6 — **DJF mean +0.9 W/m², essentially
exact.** That rules out the Pithan et al. (2014) supercooled-liquid route for midwinter.
(The Oct/Nov deficit is real and is a *separate* autumn problem.) The +30 pp polar-night
cloud excess against CERES is retrieval failure over snow at night, not evidence.

*Unresolved:* winter cover is **unchanged** (0.963 vs 0.964 — both formulas saturate at
43 cm) and winter SWE is slightly *higher*, so the DJF response is not a direct cover or
mass effect. `str` and `sshf` shifts are the surface being colder, i.e. responses.
~~`stl1` and `skt` are unusable (~0.07 K — the same dead-field trap as pressure-level `r`).~~
**That was wrong** — both are sound (210–319 K); the 0.07 was 273/3600, an analysis bug that
divided an *instantaneous* temperature by the accumulation period (fixed in `boxcache.py`).
Using them properly: G4 warms skin (+0.990) and `stl1` (+0.883) as much as T2m (+0.952), so
its warming is genuine, **but I1/I2 cool `stl4` by ~2.1 K** — and LPJ-GUESS is fed
`ST1L`–`ST4L`, so that reaches the vegetation model.
H1 cooled DJF by −1.233 through a completely different route to the same tile fractions,
so this is the second time cutting snow-tile fraction has cooled winter hard.

### Round 16 — in flight (2026-08-04)

`vlamsk_mod.F90` hardcoded the snow-tile skin conductivities, so they were untunable
until now. λ_sk sets how tightly the skin couples to the medium below; a small value lets
it radiatively decouple and crash in polar night. Physically λ ≈ k_snow/d_skin with
k_snow 0.1–0.4 W m⁻¹ K⁻¹ and d_skin 0.01–0.05 m, so **roughly 2–40**. As released:

| constant | value | tile |
|---|---:|---|
| **`ZSNOW`** | **7** | 5, **exposed snow** — the weakly coupled one, and the tile covering tundra (25.6 % of the box) |
| `ZSNOW_GLACIER` | 8 | 5, where SWE > 9000 mm |
| `ZSNOWHVEG` | 20 | 7, snow under high vegetation |
| `ZLARGESN` | 50 | 5, when melting |

| run | setting | prediction |
|---|---|---|
| **J1** | I1 + `ECE_LAMSK_SN: 15` | DJF +1…+3 K, JJA barely moves |
| **J2** | I1 + `ECE_LAMSK_SN: 25` | same, stronger |

**Falsifier on record: if DJF does not respond, skin conductivity is not the route** and
the winter bias lies in boundary-layer mixing instead (Holtslag et al. 2013, Sandu et al.
2013). Also watching that melt-out and May/June cover **stay** fixed — drift back toward
the control would mean the two changes interact and the pair is not separable.

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

**✅ RESOLVED (round 15): the link is weak, and Route B is not where the kelvin lives.**
I2 fixed the snow field to satellite accuracy — melt-out 1 day late instead of 15 — and
gained **+0.094 K** in JJA, inside the noise floor. Per W/m² of surface shortwave the snow
route is **4× less efficient** than the vegetation route (0.03 vs 0.12 K per W/m²), because
energy added over a *melting* pack goes into melt rather than sensible heat. Route B was a
real, satellite-measured radiative bias with almost no temperature leverage. **Fixing it is
still right** — the snow field is now correct, which is what matters for LPJG growing-season
onset — but it is not the remaining boreal kelvin.

**⚠ Coupled risk, now with a fix available.** G4 alone delays melt-out 24 → 28 May, which
is later growing-season onset for LPJ-GUESS and a plausible route to harming forest
establishment while AMIP T2m improves. **`ECE_SNOW_SCF=1` removes that** (melt-out 12 May
against satellite 11 May) — but on its own it costs −2.7 K in DJF until the skin-conductivity
question (round 16) is settled. Do not ship the snow fix into a coupled run until J1/J2
report.

## 3c. The near-surface cold is global, not boreal (2026-08-04)

**Already known** (report `sub:vprof`): above 700 hPa Siberia merely shares a global
tropospheric cold bias that no boreal lever can touch; the Siberia-specific excess lives
below 850 hPa. The F-series was framed correctly.

**New.** Resolving the 2 m bias by dominant surface type: **all 20 types are cold in every
season, none warm.** Global land DJF −1.73, MAM −1.60, JJA −1.35. Bare soil is the worst
(−2.16 JJA), tropical evergreen broadleaf −1.46.

Two controls bound it:
- **elevation-dependent** — −0.70 K below 200 m rising to −2.43 K above 2000 m (smoothed
  model orography vs ERA5's 25 km), so all-land figures overstate it
- **prescribed-SST ocean is −0.72 K** — ~0.7 K exists where the land surface is absent

**Extending the column into the soil** (`vertical_bias_column.py`), JJA, land <500 m:
soil L1 −0.36, L2 −0.38, L3 −0.72, L4 −0.89, skin −0.53, **2 m −1.00**, 850 hPa −0.92,
700 −0.89, 500 −1.24, 200 −2.97. **The soil and skin are the least-biased parts of the whole
column**, and 2 m sits ~0.5 K below the skin under it. The ground is not holding the
near-surface cold.

**Consequence.** Siberia's JJA −2.58 ≈ 0.7 universal + 0.3 orographic + 1.3–1.6
boreal-specific. **G4 has taken ~0.95 of the boreal budget** — which explains F5 = F4 far
better than any parameter limit. *The boreal budget is nearly spent.*

**Benign side-effect:** because every type is cold, a vegetation-indexed lever does not damage
the regions it leaks into (types 4/9 outside the boreal zone are cold too, JJA −0.68/−1.06).
That licenses the round-17 restriction.

## 3d. ⭐ NEW DIRECTION (2026-08-04): the clear-sky SURFACE shortwave deficit

With the boreal surface budget nearly spent, the dominant remaining term for land 2 m
temperature is the global tropospheric cold bias `sub:vprof` identified weeks ago and nobody
has worked on. First look (`tropo_bias_section.py`):

**Shape.** The cold **maximises at the tropopause at every latitude** — 200–300 hPa in the
extratropics, ~100 hPa in the tropics — reaching −4.2 K (90–60S) and −4.5 K (60–90N)
annually, −5.6 K over the Arctic in JJA. The lower troposphere is a fairly uniform −0.6 to
−1.2 K everywhere except 90–60S, which is near zero.

**Energetically consistent, and NOT cloud:**

| global | model | CERES | diff |
|---|---:|---:|---:|
| absorbed SW | 239.98 | 241.36 | **−1.38** |
| planetary albedo | 0.2951 | 0.2900 | +0.0051 |
| SW CRE | −44.42 | −45.32 | **+0.90** ← clouds reflect *less* |
| **clear-sky absorbed SW** | **284.40** | **286.68** | **−2.28** |

Cloud is the *wrong sign* to explain the deficit. **Every lever in 38 runs was a cloud or
surface-vegetation lever; none touched clear-sky shortwave.**

**And it is at the SURFACE, not in the atmosphere** — atmospheric clear-sky absorption is
**+0.78** (wrong sign, so not aerosol/vapour/ozone), while the surface absorbs −3.46:

| | per unit area | global contribution | Δalbedo |
|---|---:|---:|---:|
| **land** | **−6.85** | −2.00 | +0.030 |
| ocean ice-free | −1.38 | −0.89 | +0.005 |
| **sea ice** | **−8.59** | −0.57 | +0.056 |

**Caveats before anyone spends a run.** CERES *surface* fluxes are a derived product
(Kato et al. 2018) with several W/m² regional uncertainty, not direct observation. The albedo
column uses CERES clear-sky downward as a common denominator because the model does not output
clear-sky downward SW at the surface. And **this must be reconciled with `sub:albreg`**, which
found all-sky land albedo near-perfect outside high latitudes — different quantities, and one
of them may be wrong. That reconciliation is the first job, not a run.

## 4. Open problems

0. **⭐⭐ THE CURRENT BLOCKER — a −4.7 K Siberian winter bias, previously hidden.**
   Correcting the snow cover removed a compensating error: excess cover had been propping
   DJF up by ~2.7 K. Against ERA5 the control is too cold in every season (−1.95 / −1.97 /
   −2.58 / −2.69) and I1 sits at **−4.71** in DJF. It is **not** a radiation-supply problem
   (DJF downward LW is +0.9 W/m² vs CERES, essentially exact, which rules out Pithan et al.
   2014 for midwinter). With radiation right and the surface still too cold, the deficit is
   non-radiative — skin decoupling or boundary-layer mixing. **Round 16 (J1/J2) tests the
   first**; if it fails, the second (Holtslag 2013, Sandu 2013) is next.

   Related and separate: a **−12 W/m² downward-LW deficit in Oct–Nov**, unexplained.

1. **⭐ RESOLVED — HTESSEL has no sub-grid snow depletion.** Cover comes from the
   **cell-mean** depth alone, `ZCVS = min(1, d̄_cm/10)`, so it stays pinned near 1 until
   the mean itself drops below 10 cm. The model's median May cell still holds **27 cm**
   (p10 4.5, p50 27.1, p90 56.5), so the formula returns ~full cover while satellite
   shows 0.647. Real 100 km regions at 27 cm mean are patchy — ridges blown bare,
   hollows holding metre drifts. That single gap produces the late melt-out (13 days),
   June cover ~2× observed, the +0.046 June albedo and the −5.6 W/m² June SW deficit.

   **It also explains why the two previous attempts failed.** Round 13 changed the
   *threshold* 10→30 cm, which does nothing at 27–51 cm depth; it could only act where
   depths sit near the threshold (13–27 cm) — October/November, exactly where the
   response appeared. The commented-out `tanh`/`sqrt` variants share the defect
   (`tanh(27/10)=0.99`). Snow brightness is the wrong direction; energy supply is not
   short in Mar–Apr.

   **Fix implemented, opt-in:** `ECE_SNOW_SCF=1` gives
   `SCF = tanh(d / (2.5·z0·(ρ_bulk/ρ_new)^m))` (Niu & Yang 2007; Noah-MP, CLM). The
   density ratio supplies the spring depletion — aged melting snow is denser and patches
   out at a given depth — and needs **no new prognostic**. Validated offline *before*
   coding, vs Rutgers 1990–2014:

   | | satellite | as-released | z0=0.01 | z0=0.02 |
   |---|---:|---:|---:|---:|
   | Jan–Mar | 0.999 | 0.964 | 0.964 | 0.952 |
   | Apr | 0.950 | 0.964 | 0.948 | 0.893 |
   | **May** | **0.647** | 0.890 | 0.772 | **0.588** |
   | **Jun** | **0.203** | 0.380 | 0.321 | **0.180** |

   Default is `0` (as-released), so it is inert unless switched on. oifs-48r1 `28b5542`.

2. **⚠ Still open, and weaker: the March–April mass deficit.** Loss −13.7/−14.0 mm/month
   rests on ERA5 SWE, which is forest-biased; may be a second effect or partly artifact.
   And **the temperature leverage of the whole Route B is unproven** — G4 gained +0.952 K
   with snow mass unchanged and cover *worse*, so fixing cover may buy far less T2m than
   −5.6 W/m² suggests. Where it should matter is **coupled**, via LPJG growing-season onset.
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
