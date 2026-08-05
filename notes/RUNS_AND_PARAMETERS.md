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

## 0. ⭐ ROUND 17 RESULTS (2026-08-05) — K1 ADOPTED, J and L falsified

**The kept configuration is now K1 = G4 + the snow-free land albedo correction.**
Section 1 below describes G4, which K1 is built on and which remains valid.

### K1 adopted

| | Siberian JJA | DJF | net TOA 60–90N | **global T2m RMSE** |
|---|---:|---:|---:|---:|
| control | — | — | −97.737 | 1.579 |
| G4 | +0.952\* | −0.385 | −97.603 | 1.543 |
| I3 | +1.322\* | −2.439\* | **−96.017** | 1.556 |
| **K1** | **+1.036\*** (t=8.29) | −0.496 | −97.463 | **1.524** |
| K2 | +0.913\* | −0.144 | −97.525 | **1.518** |
| *CERES* | — | — | −97.98 | — |

Both clean, every season inside its own threshold, the two best global T2m RMSE of the
campaign. Tundra exclusion worked — no Arctic penalty, against I3's +1.72 W/m².
In the Siberian box neither K is significant vs G4 (+0.084 / −0.039, inside ±0.244), **as
predicted in advance**: K's mass is crops, semidesert and the subtropical desert belt.

**What it bought, on the metric it was designed for** (global land, vs G4):

| | land T2m | land SWnet | land albedo |
|---|---:|---:|---:|
| **K1** | **+0.089 K** | +0.77 W/m² | −0.0059 |
| K2 | +0.076 K | +0.48 | −0.0038 |

The lever worked: it removed **74 %** of the +0.0080 snow-free residual, at a measured
sensitivity of **0.116 K per W/m²** against the 0.12 assumed. Geography right too — land
<60N warms +0.109 K, land >60N *cools* 0.055 K.

⚠ **Prediction was too optimistic.** ~~~1.1 W/m² and 0.1–0.2 K~~ → measured 0.77 and 0.089.
The sensitivity assumption was right; the forcing estimate was not.

### 🛑 THE NEGATIVE RESULT MATTERS MORE THAN THE 0.089 K

Three quarters of a **directly measured** land albedo bias removed → 0.089 K of the ~0.7 K
universal bias. A *perfect* fix buys ~0.12 K, **one sixth** of the target.
**Land surface albedo is closed out as the explanation for the universal cold bias.**

⚠ **K3 RETRACTED before it was run.** ~~Push the soil-albedo scale further since K2 alone
gives the best global RMSE.~~ The residual 26 % is worth ~0.03 K, and 0.95 already delivers
the measured −0.0163 — going further is tuning **against** the observation, the same
overfitting-to-present-day trap that ruled out editing the soil ancillary. **Albedo tuning
stops here.**

### J falsified (round 16) — skin conductivity is not the winter route

| Δ vs control | DJF | JJA | SON | global T2m RMSE |
|---|---:|---:|---:|---:|
| I1 (base) | −2.759\* | +1.250\* | −0.766\* | 1.627 |
| **J1** (λ_sk 15) | **−2.834\*** | +1.313\* | −1.132\* | 1.581 |
| **J2** (λ_sk 25) | **−3.226\*** | +1.246\* | −0.879\* | 1.597 |
| control | — | — | — | 1.579 |

Predicted **DJF +1…+3 K**. Measured: DJF got **worse, monotonically with the knob**, and both
runs score a global RMSE worse than the control. The pre-registered falsifier fired. With
λ_sk excluded and the route known to be non-radiative, **boundary-layer mixing** (Holtslag
2013, Sandu 2013) is what remains. Doubly moot: built on I1, now deactivated.

### L falsified (round 17) — the moss proxy, monotonically the wrong way

| | JJA | vs G4 | verdict (±0.244) |
|---|---:|---:|---|
| G4 | +0.952\* | — | the base |
| L1 (λ_sk 2.9) | +0.784\* | −0.168 | inside → null |
| **L2** (λ_sk 1.7) | +0.638\* | **−0.314** | **significant degradation** |

More moss is *worse* — opposite of Gaillard (2025), and the monotonic ordering makes noise
implausible. ⚠ **Caveat:** λ_sk is only a *proxy* for an organic layer — it changes skin
coupling but not heat capacity, moisture retention or the evaporative regime. This falsifies
**the proxy**, not necessarily the mechanism. A real organic layer needs a soil-column change.

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

### ⭐⭐ I3 evaluated (2026-08-04) — strongest summer lever built, NOT adopted

I3 = G4 + mode 2 (SDOR-scaled), calibrated to match I1 **over the Siberian box** so any
remaining difference would be spatial structure. It does not behave like I1.

**What it wins:**

| | control | G4 | I1 | **I3** | CERES |
|---|---:|---:|---:|---:|---:|
| Siberian JJA T2m [°C] | 9.728 | +0.952 | +1.250 | **+1.322** (t=10.61) | — |
| Siberian surface net SW [W/m²] | 153.78 | +8.10 | +12.77 | **166.58** | **166.26** |
| global T2m RMSE vs ERA5 [K] | 1.579 | −0.036 | **+0.048** | **−0.023** | — |
| global net TOA [W/m²] | 0.643 | −0.128 | −0.021 | +0.056 | — |
| SO SW RMSE | 6.877 | −2.076 | −1.906 | −2.065 | — |

+1.322 K is the **largest JJA gain of 36 evaluated runs**, and the Siberian surface SW
deficit that opened the whole albedo thread **closes to within 0.3 W/m² of CERES**. Where I1
*degrades* global T2m RMSE, I3 improves it. Per season it dominates I1 outright:

| Δ vs control [K] | DJF | MAM | JJA | SON |
|---|---:|---:|---:|---:|
| threshold (own season) | ±0.601 | ±0.394 | ±0.244 | ±0.441 |
| G4 | −0.385 | −0.360 | +0.952\* | +0.185 |
| I1 | −2.759\* | −0.306 | +1.250\* | −0.766\* |
| **I3** | **−2.439\*** | −0.214 | **+1.322\*** | **−0.105** |

**⚠ What only the FULL guardrail set revealed** — the Siberian-box scripts alone would have
promoted I3, and the verdict "best lever of the campaign" *was* drafted on them before the
guardrails were run:

| guardrail | control | I1 | **I3** | CERES |
|---|---:|---:|---:|---:|
| net TOA 60–90N [W/m²] | −97.74 | +0.87 | **−96.02 (+1.72)** | **−97.98** |
| net TOA NH [W/m²] | 0.280 | +0.245 | **0.952 (+0.672)** | **0.08** |
| NH−SH albedo | −0.003 | −0.002 | **−0.007 (−0.004)** | — |
| Nordic Seas SW RMSE | 9.058 | +0.196 | **9.557 (+0.499)** | — |

This is the mechanism **working as designed** — less high-latitude snow cover admits more SW
at high latitude — but the Arctic had no energy spare. The control sat almost exactly on the
observation at 60–90N (−97.74 vs −97.98). The NH TOA and NH−SH albedo changes are each the
**largest of all 37 runs**, and Nordic Seas is the **only significant SW RMSE degradation in
the campaign** (`rmse_significance.py` flags two runs there; I3 is one, +0.07 vs ±0.051).

*Caveat, so this is not over-read:* these are PI runs scored against present-day CERES, and
the standing rule is that the PI target is net TOA ≈ 0, **not** the CERES column. That
weakens the "away from CERES" reading of the two regional TOA rows. It does **not** weaken
NH−SH albedo or Nordic Seas, which are model-vs-model and model-vs-RMSE.

**Why this matters more for coupled than AMIP:** all four are Arctic energy-*gain*
guardrails, and AMIP structurally hides the consequence — prescribed SST and sea ice absorb
1.7 W/m² without responding. The coupled PI spin-up will not.

### ✅ The Siberian winter penalty is NOT radiative (2026-08-04)

| Siberian DJF | snow albedo | SWE [m] | absorbed SW [W/m²] |
|---|---:|---:|---:|
| control | 0.6789 | 0.6933 | 8.16 |
| G4 | 0.6760 | 0.6901 | 8.18 |
| I1 | 0.6800 | 0.6926 | 8.16 |
| I2 | 0.6796 | 0.6958 | 8.16 |
| **I3** | 0.6797 | 0.6944 | 8.22 |

Albedo, snowpack and absorbed SW are **unchanged to three decimals**, on a midwinter budget
of only ~8 W/m². **A −2.4 K penalty cannot be built from that.** This generalises the finding
that killed the compensating-error story (winter cover 0.963 vs 0.964) from cover to the
whole radiative chain, and it holds for I2, which carries the penalty with no vegetation
change at all. The route is **thermal — the snow/exposed-tile coupling** — which is exactly
what round 16 tests.

### 🛑 `ECE_SNOW_SCF` DEACTIVATED (2026-08-05) — and what reactivation needs

**The switch stays at 0 (as-released). No production run sets it.** The default is
`ECE_SNOW_SCF: 0`, so this is a decision not to use an option, not a code removal — the
implementation stays in the source for whoever picks the problem up.

**Why — direct observations, not ERA5.** RIHMI-WDC (Sherstyukov v3) station soil temperature,
43 stations in the box, **under natural cover** so the snowpack is intact. Sampled at the
station points, not box-averaged:

| | obs depth | OBS | model | bias | OBS offset | model offset |
|---|---|---:|---:|---:|---:|---:|
| `stl2` | 20 cm | −6.07 | −5.03 | **+1.04** | **+17.78** | +20.11 |
| `stl3` | 40/80 cm | −3.58 | −3.19 | **+0.39** | **+20.28** | +21.95 |
| `stl4` | 120–240 cm | +0.34 | −0.23 | **−0.57** | **+24.18** | +24.85 |

Offset = soil − 2 m air, each dataset using its **own** air (insensitive to the 2.8 K air
change the I-series also carries). Scored on offset:

| | stl1 | stl2 | stl3 | stl4 |
|---|---:|---:|---:|---:|
| control | 19.87 | 20.39 | 22.24 | 25.13 |
| G4 | 20.27 | 20.81 | 22.72 | 25.82 |
| **I1** | **1.53** | **2.61** | **6.36** | **15.25** |
| **I3** | 1.41 | 2.54 | 6.46 | 15.63 |
| **OBSERVED** | — | **+17.78** | **+20.28** | **+24.18** |

I1 is **−15.2 K wrong at 20 cm**. Raw: its soil sits at −25.8 °C where stations measure
−6.1 °C. **The control was correct to ~1 K at every verifiable depth**, so this is a NEW
error, not the removal of an old one. `ST1L`–`ST4L` feed LPJ-GUESS → disqualifying.

⚠ **`stl1` (0–7 cm) has NO winter observations** — the Russian network withdraws shallow
Savinov thermometers for the cold season. `stl2` is the shallowest verifiable layer.

⚠ **Correction:** the "+10 to +20 K" offset band quoted earlier was from memory of the
permafrost literature and is too low. The station record gives **+17.8 to +24.2 K**.

**The mechanism is UNIDENTIFIED**, and the obvious candidates are ruled out:
- ~~the tanh cuts cover in October when snow is shallow~~ — **wrong**: cover changes −0.02 in
  Oct, SWE not at all.
- The soil **leads** the air (Oct: soil −4.19 K, air −0.02) → surface coupling, not atmosphere.
- Winter snow density falls 20–27 kg/m³ with unchanged SWE → a *deeper, less dense* pack,
  which insulates **better**. Density is a consequence of a cold pack, not its cause.

**Reactivation needs three things, in order:**
1. **Diagnose `ZCVS`** — it is not in the output stream, which is exactly why the mechanism is
   still open. Rebuild nothing before the snow-tile fraction can be seen.
2. **Accumulation/ablation hysteresis** — we apply one curve year-round; Niu & Yang (2007) and
   Swenson & Lawrence (2012) both use different curves for accumulation and ablation. Matches
   the measured timing: benefit in May–June, damage from October, cleanly separable in season.
3. **Re-verify on the SOIL OFFSET, not snow cover.** Round 15 accepted it because melt-out
   matched satellite to within a day. That was true and insufficient.

**Weigh the prize first:** I2 (scheme alone) bought **+0.094 K JJA — inside the noise floor**.
The value is a correct snow field for LPJG phenology, not temperature. That does not justify a
source change plus a run family while a 20 K soil error is unexplained.

*Asset:* `soil_temp_vs_rihmi.py`; obs at `/work/ab0246/a270092/obs/RIHMI-WDC/`.

### Round 19 — N series, DAILY snow/soil process diagnostic, IN FLIGHT (2026-08-05)

⚠ **NOT tuning candidates.** N1/N2 exist to test two inferences that monthly data cannot
settle — and that I had to **retract once each** on 2026-08-05 before building on them.
K1 base, 10 yr (1870–1880), identical except the scheme under test.

| run | change | asks |
|---|---|---|
| **N1** `amip_N1_snowdiag` | K1 + daily `sd`/`rsn`/`tsn`/`asn`/`stl1`/`stl2` | reference |
| **N2** `amip_N2_snowdiag_scf` | N1 + `ECE_SNOW_SCF: 1`, `Z0: 0.016` | the scheme, at daily resolution |

**INFERENCE 1 — the winter soil collapse is SEEDED IN AUTUMN.** `ZCVS` is not an output
field (it equals `FRTI(5)+FRTI(7)`, and neither is written either), so it was
reconstructed offline from each run's own snowpack:

| month | Sep | **Oct** | Nov | Dec | Jan |
|---|---:|---:|---:|---:|---:|
| ΔZCVS | +0.026 | **−0.075** | −0.008 | −0.001 | −0.000 |
| Δsoil [K] | +0.0 | **−4.2** | −16.2 | −22.2 | −23.8 |

The cover difference **peaks in October exactly as the soil starts to diverge**, then
vanishes while the soil runs to −24 K. Reading: ~7 % more bare ground in October, before
the insulating pack closes, cools the soil while it still can be; from November the pack
seals and the anomaly is locked in. Midwinter cover is **identical** (Jan ΔZCVS −0.0004).
⚠ **NOT DEMONSTRATED**: that the Nov amplification (−4 → −16 K) follows from the Oct seed.

**INFERENCE 2 — would a melt-state gate fire when the spring depletion is needed?**
Monthly-mean `tsn` is 271.0 K in May, reaching 273.0 only in June, yet melt-out is late
May. If the pack is ripe only intermittently in May the gate fires **late**, the spring
depletion is lost, and the scheme becomes all winter cost and no summer benefit.

**WHY THE HYSTERESIS WAS NOT SIMPLY BUILT** — both found by trying to falsify it:

1. **A ripe/not-ripe STATE TEST IS NOT HYSTERESIS.** It is reversible: cover would flip
   back to the linear formula every night and every cold snap, giving a spurious diurnal
   oscillation in snow-covered fraction. Real melt-out is irreversible within a season.
   The correct formulation (Swenson & Lawrence 2012) carries a prognostic seasonal
   `SWE_max` and depletes against it — but that is a **new surface prognostic**.
2. ~~**Density cannot substitute for melt state.** October snow is ρ≈170, May ρ≈315, but
   May depth is still 0.29 m — so any exponent that saturates cover in October also
   saturates it in May. The two regimes are inseparable by density at their respective
   depths.~~ **WRONG — corrected the same day, see the O-series entry below.** That was
   argued from the *current* parameters, not from what the formula can do. Density
   separates them fine (m>1.23 suffices for the right ordering, m=4 gives a factor 5.5);
   m=1.6 merely does it so weakly (separation 0.077) that the October cover deficit is
   the price paid for a negligible spring gain. This is what the O series tests.

**WHAT TO LOOK FOR**
- *Sep–Nov*: does `stl1` diverge only on days the cover differs, or keep falling after the
  cover difference closes? Only the latter confirms seeding + memory.
- *Apr–Jun*: fraction of days with `tsn` ≥ 273.15. If small in May, inference 2 fails and
  ripeness is the wrong criterion.

**IMPLEMENTATION NOTE.** Daily output added via a **per-run** `add_config_sources`
override (`namelists/oifs/48r1/xios/file_def_snowdiag.xml.j2`), so the shared template is
untouched and no other run is affected.

⚠ **Two process traps hit while building this**, both worth avoiding again:
- The first patch landed in a **6 h** block: the remapped section has one *before* the
  daily block, and "first anchor after `atm_remapped`" found the wrong one.
- The first verification **counted `field_ref="sd"` across the whole file** — already 12
  from pre-existing `1mo`/`6h` blocks — so it reported success on a no-op. **Verify
  per-block, in the GENERATED `file_def.xml`, not the template.**

### Round 19 — O series, the RE-PARAMETERISED depletion, IN FLIGHT (2026-08-05)

Paired with N1/N2 above: same K1 base, same 10 yr, same daily output, so all four
difference directly. **Namelist-only** — no source change, no new prognostic.

| run | z0 | ρ_new | m | predicted SCF Oct | predicted SCF May |
|---|---:|---:|---:|---:|---:|
| N2 (current) | 0.016 | 100 | 1.6 | 0.894 | 0.817 |
| **O1** `amip_O1_scf_m4` | 0.018 | 170 | **4.0** | **0.995** | **0.495** |
| **O2** `amip_O2_scf_m3` | 0.018 | 170 | **3.0** | 0.995 | 0.764 |

**⚠ CORRECTION to the N-series entry above.** It records that "density cannot substitute
for melt state ... the two regimes are genuinely inseparable by density at their
respective depths". **That is wrong**, and was based on the *current* parameters rather
than on what the formula can do. Density separates them fine; m=1.6 simply barely does.

**The arithmetic.** SCF = tanh(d/L), L = 2.5·z0·(ρ/ρ_new)^m. The ratio of the two tanh
arguments for the two months that matter, using Siberian box means (d_Oct 0.135 m,
ρ_Oct 170; d_May 0.288 m, ρ_May 315):

```
x_Oct / x_May = (d_Oct/d_May) · (ρ_May/ρ_Oct)^m = 0.469 · 1.853^m
```

**Depth works against us** — October's pack is half as deep, pushing it toward *less*
cover. Density works for us, and whether it wins is entirely down to `m`:

- **m > 1.23** needed merely for October to hold more cover than May
- **m = 1.6** clears that bar but gives a ratio of only **1.23** → SCF 0.894 vs 0.817,
  separation **0.077**. Almost none.
- **m = 4** gives **5.5** — ample.

So the fault is not that density cannot discriminate; it is that at m=1.6 it *barely*
does, and the **−0.075 October cover deficit that seeds the winter soil collapse is the
price paid for a spring depletion only 0.08 deeper**.

**Raising m alone makes October worse.** October ρ=170 already exceeds ρ_new=100, so it
sits on the depleting branch and a steeper exponent depletes it harder (m=4, ρ_new=100 →
SCF_Oct **0.383**). Both parameters must move: recentre ρ_new on autumn density so
October becomes the reference state, *then* steepen.

**Why this beats the hysteresis it replaces.** It sidesteps the reversibility objection
that killed the ripe/not-ripe state test: **density is a genuine physical memory
variable**. It ripens monotonically through the season and does not flip back overnight,
so the depletion cannot walk backwards the way a melt-state switch would. And it needs
no new surface prognostic.

**⚠ This is curve-fitting until the daily data supports it:**
- m=4 is far from Niu & Yang's calibrated 1.6 and needs justifying as more than a fit to
  the two numbers it was derived from.
- Those are **box-mean monthly** densities; the model sees a cell-by-cell distribution and
  a steep exponent **amplifies spread**, so the box mean may mislead badly. N1/N2's daily
  output is what settles this.
- Recentring ρ_new changes depletion **globally**, including where the original
  calibration was validated against Rutgers/IMS satellite cover.

**PRE-REGISTERED FALSIFIERS**
1. **Winter soil must return to the N1 reference.** If O1 still shows the −16 K DJF
   collapse, the density route is dead and the prognostic `SWE_max` rewrite is the only
   remaining option.
2. **May/June cover must still deplete.** If not, the summer gain is lost and O1 is all
   cost — the same failure mode predicted for the ripeness gate.

### Round 18 — aerosol diagnostics, IN FLIGHT (2026-08-05)

⚠ **On the PRESENT-DAY base (1989-2015), not the PI base.** MACv2-SP anthropogenic
aerosol is TRANSIENT, so at the usual 1872-1915 the plumes are already near zero and
these tests would show nothing. Present-day is also the period the +2.68 W/m² was
measured over and the only one comparable to CERES. **Listed in `NOT_LEVERS`** so
`evaluate.sh` cannot difference them against the PI control over the wrong years.
Compare against `amip_presentday`.

| run | change | asks |
|---|---|---|
| **M1** `amip_M1_noanthaer` | `LMACV2SP: false`, `LMACV2SP_CCNF: false` | is the clear-sky excess ANTHROPOGENIC? |
| **M2** `amip_M2_aer3d` | `LAER3D: true` | is it the VERTICAL DISTRIBUTION of the same aerosol? |

**The aerosol configuration, traced end to end** (several wrong turns before this stuck):

| component | source | year-dependent? | tunable? |
|---|---|---|---|
| anthropogenic | MACv2-SP plumes, `SPv2.1_1850-2023_CMIP7.nc` | **yes** | plume amplitude |
| natural background | CAMS/MACC monthly climatology (`NAERMACC=1`) | **no, fixed** | **dust only** |
| GHG, ozone, strat. aerosol | CMIP7 | yes | — |
| `CMIP6_..._aerosol_radiation_2D` file in the input dir | **not read by this version** | — | — |

- `suecrad.F90:1115`: `IF (LMACV2SP) NAERMACC = 1` — CAMS supplies the background,
  MACv2-SP adds the anthropogenic plumes on top.
- `LAER3D=.false.` comes from `esm_tools/configs/components/oifs/oifs48.tuning.yaml`.
  **ECMWF's own TCO95L91 reference `fort.4` sets nothing in `&NAERAD`**, so ECMWF takes
  the shipped default `.TRUE.` (3D CAMS). Ours is an EC-Earth/AWI choice.
- The 3D file **is** staged (`ifsdata/aerosol_cams_climatology_43R3_3D.nc`) and carries
  the SAME mass as the 2D to 0.4 % — so M2 changes only *where* the aerosol sits.
- **Sea salt has NO scaling knob anywhere in the radiation tree.** Only dust does
  (`RDUMULF`, `RWGHTDU1/2/3`), and those are **not namelist-settable** — hardcoded in
  `suecrad.F90` (`LAERADJDU=.FALSE.`, all factors 1.0), so no aerosol scaling is active
  in any campaign run.
- ⚠ `LMACV2SP_CCNF=.true.` feeds aerosol into cloud droplet number, so **neither run is
  a clean clear-sky lever** — score both columns.

**Expectations, and the falsifier.** MACv2-SP is NH continental outflow, not a Southern
Ocean signal. If M1 recovers most of the +4.22 at 30-60N but leaves the SO's +5.10
intact, the residual is the CAMS sea-salt background — and `Sea_Salt_bin2` is 77 % of
SO AOD (119 mg/m², implied AOD550 0.107 of 0.140).

**Prior evidence says expect little.** MISR gives SO AOD 0.107 vs the model's ~0.140;
that +0.033 is worth only ~0.7 W/m², **13 % of the +5.10**. And the clear-sky deficit
is predominantly a SURFACE problem: globally the surface absorbs 3.46 W/m² too little
while the atmosphere absorbs 0.78 too *much*.

⚠ **DO NOT ADOPT EITHER.** Removing anthropogenic aerosol from a present-day run is
physically wrong; these attribute a residual, nothing more.

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

⚠ **Both J runs are built on I1, before I3 was shown to be the better base in every season.**
That was correct with the information available, and the round is still diagnostic — it tests
whether λ_sk is the winter route at all, and I1/I3 share the same non-radiative penalty, so
an answer on I1 transfers. **But if λ_sk proves to be the route, the production pair must be
re-run on I3**, and the Arctic guardrails re-checked on the combination rather than assumed
to superpose. Superposition has been wrong in *sign* twice in this campaign.

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

## 3e. ⭐ Land albedo (2026-08-04): half of it is the SNOW TILE

*This is the reconciliation §3d asked for, and it mostly vindicates `sub:albreg`.*

Measured assumption-free (each dataset divided by its **own** downward flux — an earlier pass
that borrowed the CERES denominator inflated the land excess from +0.015 to +0.030), all-sky
land surface albedo is **+0.0154** too high.

**The trap: snow is a TILE, not a vegetation type.** The first attribution split that by
dominant vegetation type from `tvh`/`cvh`/`tvl`/`cvl`. That map contains no snow — in HTESSEL
snow is a separate tile (5 on low veg, 7 under high) that **overlies** every vegetation type.
So each type was part-time snow and the ranking was largely a map of *snow duration*: evergreen
needleleaf 53.9 % of the year, tundra 63.2 %, "bare soil" 49.2 % (Sahara and high-Arctic bare
ground in one bucket). An annual-mean T2m > 5 °C filter does **not** fix this — a cell
averaging 6 °C still has months of snow.

Masking snow **per cell per month** (SWE < 1 mm, sunlit months only) splits it almost exactly
in half: **snow +0.0074, snow-free surface +0.0080.** The per-type ranking inverts:

| type | snow-covered | all months | **snow-free** |
|---|---:|---:|---:|
| Evgr Needleleaf | 53.9 % | +0.0250 | **+0.0032** |
| Bogs/Marshes | 37.1 % | +0.0143 | +0.0039 |
| Evgr Shrubs | 41.2 % | +0.0150 | +0.0058 |
| Evgr Broadleaf (tropical) | 6.9 % | +0.0017 | **−0.0003** |
| **Crops** | 25.9 % | +0.0215 | **+0.0204** |
| **Semidesert** | 22.5 % | +0.0342 | **+0.0169** |
| **bare soil** | 49.2 % | +0.0240 | **+0.0145** |
| **Irrig Crops** | 24.6 % | +0.0253 | **+0.0138** |
| **Tundra** | 63.2 % | +0.0064 | **+0.0105** ← gets *worse* |

**Every forest type falls inside +0.006 once snow is removed**, tropical broadleaf to −0.0003.
So `RVVEGALB`'s high-vegetation entries are **sound** and `sub:albreg`/`sub:albsnow` hold. What
is *new* is the +0.0080 residual on **sparse and cultivated** surfaces where the soil background
shows through — plus tundra, where snow was hiding a genuine snow-free *summer* error.
**Nothing in 41 runs has touched that half.**

**The snow half already has a lever** (model-internal, snow-covered sunlit months):

| run | snow-month albedo | vs G4 | all-land |
|---|---:|---:|---:|
| G4 | 0.5185 | — | 0.2586 |
| I1 = G4 + mode 1 | 0.5011 | −0.0174 | −0.0024 |
| I2 = mode 1 alone | 0.4996 | −0.0189 | −0.0018 |
| **I3 = G4 + mode 2** | **0.4689** | **−0.0496** | **−0.0052** |

**I3 removes ~70 % of the snow half**, I1 ~32 % — and I3 was calibrated to match I1 *in the
Siberian box*, so globally it is **twice as strong as intended**. Same over-reach the Arctic
guardrails bill it for: the SDOR scaling does most of its work outside the box it was tuned in.

### ⚠ The two halves cannot simply be stacked

Keeping I3 and adding a vegetation-albedo correction on top is **not additive in the direction
that matters**. A darker land surface at high latitude adds absorbed SW to exactly the region
I3 has already pushed 1.72 W/m² past CERES. They compete for one energy budget.

- **Mid-latitude / subtropical entries are the safe ones** — crops (+0.0204), semidesert
  (+0.0169), bare soil (+0.0145), irrigated crops (+0.0138) lie outside the Arctic, so fixing
  them attacks the ~0.7 K **universal** bias without adding to the 60–90N surplus.
- **Tundra (+0.0105) is the dangerous one** — the single entry that *would* stack with I3,
  being boreal and in-season, hence both the most attractive and the most likely to break the
  Arctic guardrails.
- **Where bare fraction dominates, `RVVEGALB` cannot reach it** — the background soil albedo
  climatology carries it, and that is an ancillary field, not a namelist parameter.

*Asset:* `land_albedo_snow_split.py`, chained into `evaluate.sh --obs`.

## 4. Open problems

0. **⭐⭐ THE CURRENT BLOCKER — a −4.7 K Siberian winter bias accompanying the snow fix.**
   ~~Correcting the snow cover removed a compensating error: excess cover had been propping
   DJF up by ~2.7 K.~~ **Retracted** — winter cover is unchanged (0.963 vs 0.964); see §3.
   Now **positively excluded as radiative at all**: DJF snow albedo, SWE and absorbed SW are
   identical across control/G4/I1/I2/I3 to three decimals, on an ~8 W/m² midwinter budget
   that cannot build −2.4 K. The route is **thermal**.
   Against ERA5 the control is too cold in every season (−1.95 / −1.97 /
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

---

## 6. ⭐ THE BOREAL FOREST FAILS BY TWO ROUTES (2026-08-05)

Evaluated from `Tuning_test_07*`/`08*` in `/work/bb1469/a270270/runtime/awiesm3-v3.4`
(`lai.out`, `fpc.out`, `est_limits.out`). Larch zone 60–72N × 100–160E, 2520 cells:

| | 080a CRUNCEP base | 09A newSeaIce |
|---|---:|---:|
| `mTmin20` | −38.94 °C | −41.89 °C |
| **BNE** gate `tcmin_est` −30 | **0.7 % clear** | **0.0 % clear** |
| **BNS** cold gate | **none → 100 %** | **none → 100 %** |
| `agdd5` vs BNS `gdd5min_est` 350 | **83.6 % clear** | 79.0 % |
| BNE LAI | **0.000** | **0.000** |
| BNS LAI | 0.129 | 0.103 |
| C3G LAI | 0.751 | 0.510 |
| **grass : larch** | **5.8×** | 5.0× |
| cells with C3G > BNS | 65 % | 60 % |

- **BNE/BINE = climate gate.** 0.7 % of cells clear −30 °C; LAI exactly 0.000. The
  winter-cold story applies **here** — evergreen west Eurasia and N America.
- **BNS (larch) = no gate, outcompeted.** No cold limit, clears GDD on 84 % of cells, so it is
  **climatically permitted** — yet grass carries 5–6× its LAI and larch is absent from half
  the cells where it is allowed. **NE Siberia is ~80 % larch, so the winter-gate mechanism
  does not explain the region that matters most.**

### Reconciliation with the falsified "competition, not climate"

| claim | status |
|---|---|
| competition is an LPJG **parameter error** | **dead** (forcing-transfer reproduces the loss unchanged) |
| competition is the **pathway** by which climate removes larch | **supported** (permitted on 84 %, loses 5.8:1) |

Forcing is the **driver**, grass-vs-larch the **mechanism**. Both hold. The report previously
recorded only the falsification, so a reader would wrongly conclude growing-season competition
had been ruled out — it had never been tested.

No conflict with the "cold-season not growing-season" seasonality result either: that measures
vegetation→climate, this measures climate→vegetation. Opposite directions of one coupling.

⚠ **Knob evidence is weak, flagged as such.** 08G (`C3G pstemp_low`=12) moves BNS 0.129→0.156
and TREEFPC 0.060→0.072 while C3G falls 0.751→0.583 — right direction, but the 08 runs sit on
different atmospheric bases (06T vs 06V), so it is not a clean attribution. On the old forcing
**07A (BNS `greff_min` relaxed) did not recover BNS at all** (0.024→0.021). A
forcing-consistent re-spin remains the precondition for competition tuning.

---

## 7. PARKED (2026-08-05) — sea ice, melt ponds, and the clear-sky SW split

**Not a lever for either campaign target.** Filed because it is a genuine model defect
that matters for sea-ice realism in the coupled PI spin-up, and because the diagnosis
cost real effort. Deliberately kept out of the tex.

### Why it is not usable for Siberia or the Southern Ocean

- **Siberia: no.** The Siberian target lives in AMIP, where sea ice is *prescribed* and
  FESIM does not exist. Coupled, the chain ice albedo -> Arctic Ocean energy -> advection
  into Siberia is long and weak against a land-surface bias.
- **Southern Ocean: marginal.** Reducing Antarctic ponds brightens the ice, which is the
  *correct* direction (surface −0.032 too dark, ponds 2–3× too abundant — consistent),
  and more reflection helps the SO energy target. But it is worth ~0.1–0.2 W/m² globally
  against a **−7.65 W/m² SO cloud deficit** that dominates the band.

### The solid finding: melt pond fraction is 2–3× too high in BOTH hemispheres

Measured from `apnd` in the gen-10 coupled runs (`/work/bb1469/a270270/runtime/awiesm3-v3.4`):

| | model | observed | ratio |
|---|---:|---:|---:|
| Arctic July | **0.51** | 0.15–0.25 (MODIS, Rösel et al. 2012) | 2–3× |
| Arctic August | 0.53 | declining | — |
| Antarctic February | **0.12** | <0.05 | 2–3× |

**The global scheme discriminates the hemispheres correctly on physics alone** — SH ponding
is a quarter of NH, and SH ice keeps 0.28 m of snow through summer while NH snow goes to
~0.001 m. **No hemispheric parameter is needed or justified**: inventing one would fix the
model to present-day geography and destroy its ability to respond in another climate state
(the same principle that ruled out a spatial land-albedo correction, §3e).

**Prime suspect: `rfracmax = 1.0`** in `&meltpond` — 100 % of meltwater retained in ponds at
maximum, against CICE's default **0.85**. Second candidate `pndaspect = 0.8`. Both global,
both physical, and a single fix should improve both hemispheres at once.

### The sea-ice albedo bias is ARCTIC, and the hemispheres oppose

AMIP (`amip_presentday`, OIFS `RALBICE_AR`/`_AN` tables — no FESIM):

| | season | model | CERES | bias |
|---|---|---:|---:|---:|
| **Arctic** | melt (JJA) | 0.4948 | 0.3752 | **+0.1196** |
| Arctic | annual | 0.5839 | 0.4810 | +0.1029 |
| **Antarctic** | melt (DJF) | 0.2490 | 0.3263 | **−0.0772** |
| Antarctic | annual | 0.4088 | 0.4412 | −0.0324 |

⚠ **Correction:** the "+0.0327 sea-ice albedo excess" quoted earlier was a *global* sea-ice
mean that averaged a large Arctic positive against an Antarctic negative and hid both.

⚠ This AMIP bias is an **OIFS table problem with no FESIM in the configuration**. It does
*not* constrain the pond scheme, and the pond scheme does not explain it.

### Two corrections to earlier claims in this thread

1. ~~Coupled runs have an output gap: only `ci` and `ssrd` are written, so sea-ice albedo is
   unmeasurable.~~ **Wrong.** Full SW output exists under CMIP names in the `atmos_*` streams
   (`rsus`, `rsds`, `rsuscs`, `rsdscs`, `rsdt`, `rsut`, `rsutcs`). Only the
   `atm_remapped_1d_*` stream is sparse. No esm_tools change is needed.
2. Coupled sea-ice albedo *is* measurable (`rsus/rsds`), but the comparison is **confounded
   by ice extent**: each run samples a different cell set, so the CERES value itself shifts
   between runs (NH JJA 0.327 for 10A vs 0.334 for 10B). A common ice mask and period are
   required before quoting any coupled albedo bias.

### The clear-sky / cloud split that framed all of this

TOA reflected SW vs CERES, model − CERES [W/m²], + = model reflects more:

| band | clear-sky | cloud |
|---|---:|---:|
| 60–90N | +4.42 | +1.36 |
| 30–60N | +4.22 | −4.22 |
| tropics | +1.55 | +0.41 |
| **SO 45–65S** | **+5.10** | **−7.65** |
| 60–90S | +2.34 | −8.19 |
| GLOBAL | +2.68 | −1.17 |

**Two independent errors that hide each other**: in the SO they nearly cancel (+5.10 and
−7.65 give only −2.54 all-sky), so the SO looks nearly right in the all-sky total while both
halves are badly wrong. Probably why the campaign's two targets have fought each other —
A1a bought 82 % of the SO cloud gap and cost −0.75 K in Siberia.

SO 45–65S clear-sky excess by surface: ice-free ocean +3.56 (79 % of area), **sea ice +11.22
(19.8 % of area, 44 % of the total)**, land +5.93 (1.2 %).

Attribution of the SO +5.10, as far as it goes:

| term | W/m² | share | evidence |
|---|---:|---:|---|
| sea-ice albedo | ~2.2 | 44 % | but see the hemispheric split above |
| ocean surface albedo | ~0.9 | 18 % | +0.0072 to +0.0109 at 45–65S, and SH-specific |
| aerosol (+0.033 AOD vs MISR) | ~0.7 | 13 % | MISR 0.107 vs model ~0.140 |
| **unexplained** | **~1.3** | **25 %** | — |

⚠ The ocean-albedo term is **asymmetric** (+0.0109 at 55–65S vs +0.0027 at 55–65N), so a
pure zenith-angle formula error is excluded — something SH-specific (whitecaps under
stronger winds? marginal-ice leakage? CERES retrieval) is involved and is not identified.

⚠ The Antarctic shows +4.61 W/m² clear-sky excess while its surface is too *dark*. That
cannot be surface albedo and is unexplained.

*Assets:* MISR/MODIS AOD are on `/pool/data/ICDC/atmosphere/{misr,modis_*}_aerosol` — no
download needed. CAMS aerosol climatology is `ifsdata/aerosol_cams_climatology_43R3*.nc`;
the 2D and 3D versions carry the same mass to 0.4 %.
