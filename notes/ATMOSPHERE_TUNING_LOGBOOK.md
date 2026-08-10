# Atmosphere tuning logbook — OpenIFS 48r1 (AMIP)

**Two things in one file:** §0 is the running record of what has been tried and what it
did; §1–10 are the reference — every knob, its default, what it does, and which are dead
ends. Started 2026-07-28 from three parallel sweeps (OIFS source, SMHI GitLab EC-Earth
repos, published literature). Add each experiment's outcome to §0 as it lands.

---

## 0. Experiment log

All AMIP TCO95, 1850 GHG, observed SST 1870s, 6 yr, evaluated **1872–1875** against
control `amip_pi_base`. Runs live in `/work/bb1469/a270092/runtime/oifsamip-cy48/`.
Evaluation script: `scripts/analysis/eval_round10_A.py`.

### Targets at the start of round 10
| | control | observed (CERES) | gap |
|---|---:|---:|---:|
| global net TOA | +0.533 | ~0 expected for PI | **+0.53** |
| global surface flux | +0.383 | 0 | **+0.38** |
| SO 45–65S TOA SW CRE | −60.06 | −68.14 | **+8.08** |
| SO cloud area | 83.15 % | 89.72 % | **−6.6 pp** |
| Siberia JJA surface net SW | 151.9 | 166.3 | **−14.4** |
| Siberia JJA cloud area | 79.1 % | 69.6 % | **+9.5 pp** |
| Siberia JJA T2m bias | — | — | ~~≈ −2.2 K~~ ~~≈ −1.0 K~~ **−1.3 to −2.0 K** |

*Corrected twice. On 2026-07-29 the −2.2 K was cut to ~−1.0 K on the argument that ~1.1 K
was reference-period offset, estimated indirectly from HadCRUT5. On **2026-07-30** the
`amip_presentday` run measured that offset directly and found only **+0.42 K**, so the
boreal target is **−1.3 to −2.0 K** — much closer to the original figure. The SO SW CRE
row is likewise **restored to ~+7.8 W/m²**: its measured period offset is −0.07, not −1.43.
See §"Period offset, measured directly" below, which supersedes the estimate.*

### ~~⚠ The target itself — about half of the boreal "bias" is the reference period~~
### Period offset, measured directly (2026-07-30) — supersedes the estimate above

`amip_presentday` is the same model and configuration run over 1989–2015, so the epoch
offset can be **measured** instead of chained through observations.

| arm | Sib JJA T2m | vs ERA5 90–14 | SO SW CRE | vs CERES | net TOA |
|---|---:|---:|---:|---:|---:|
| control (1872–1915) | 9.73 | **−2.45** | −60.29 | **+7.85** | 0.64 |
| presentday (1990–2014) | 10.14 | **−2.04** | −60.36 | **+7.78** | 2.20 |
| **measured offset** | | **+0.42 K** | | **−0.07 W/m²** | +1.56 |

**Both of yesterday's corrections were too large.** The boreal offset is **+0.42 K**, not
+1.12 K, so the period-clean bias is **−2.04 K** — close to the original figure, *not* ~1.0 K.
The SO SW CRE offset is **−0.07**, not −1.43, so that gap is **+7.78**, essentially the
original +8.08. The "~6 W/m²" claim is withdrawn.

**Why the estimate failed, and the real finding underneath.** It applied an **observed**
epoch change (~+1.1 K from HadCRUT5×amplification) to a model that does not reproduce it.
The model warms Siberian JJA by only +0.42 K between the epochs, so **the model
under-warms the historical period by ≈0.7 K** — a *trend* error distinct from the mean-state
bias, and a new result in its own right.

**Target.** If observations are right about the epoch change, the 1870s bias is ≈**−1.3 K**
(12.18 − 1.12 = 11.06 target vs 9.73). By the model's own internal offset it is ≈**−2.0 K**.
So the honest range is **−1.3 to −2.0 K**.

*Caveat:* `presentday` global net TOA is **+2.20** against CERES's +0.97 — a +1.23 excess,
larger than the control's +0.64. Aerosols were verified correct (MACv2-SP transient, years
1989→2002 logged), so this is not a missing-forcing artefact but a genuine radiative bias,
consistent with the SO cloud deficit being identical across epochs (+7.85 vs +7.78).

### Run inventory — every directory under `runtime/oifsamip-cy48/`, including the dead ones

There are **42** experiment directories but only **30** carry data. Recording which are which,
so nobody has to guess whether an empty directory is a failure worth investigating.
*(Counts change as rounds are added — regenerate with the one-liner at the end of this section.)*

**Carrying evaluated data (30).** control `amip_pi_base` (50 yr) and `amip_picontrol` (48);
`amip_presentday` (26); the levers `A1a A1b A1c A2_kknumland150 expA B1–B8 AB ABB8 C1 C2 E1`
(46 yr each); round 11 `D1 D2a D2b` (48 yr each); round 12 **`F1–F5` (48 yr each, COMPLETE)**.
All appear in the results tables above.

**Running (2).** Round 13, the H-series, testing the snow-cover-fraction lever `RQSNCR`
(critical snow depth 10 cm → 30 cm): `amip_H1_snowcr30` (the change alone, vs control) and
`amip_H2_G1_snowcr30` (the change on top of G1). Started 2026-08-02.

**Completed since this inventory was written (1).** `amip_G1_F4_D2b` — now in `RUNS` and in
all evaluator tables; see the G1 section above. Adding a run means editing
`scripts/analysis/runs.py` only, not three separate copies.

**Superseded, no output (1).** `amip_A2_kkland150` — a first A2 attempt from 2026-07-29,
replaced by `amip_A2_kknumland150`. Two runscripts exist for A2 (`_A2_kkland150.yaml` and
`_A2_kknumland150.yaml`); **only the `kknumland150` one was ever evaluated.**

**LPJG forcing generation, not tuning runs (11).** `amip_lpjgforce_chk`, `amip_nolpjg_forc`,
`amip_nolpjg_forcing`, `amip_nolpjg_pi1870`, `amip_pi_clean1/2`, `amip_pi_dbg1/2/3`,
`amip_pi_fixtest`, `amip_pi_forcing` — all from `_lpjgforcing.yaml`, 3–5 July 2026, during
the work on the clean AMIP–noLPJG forcing generator (report §`sub:amipgen`). They were built
to emit the **daily LPJG forcing set**, not the standard evaluation fields, which is why they
have zero `atm_remapped_1m_2t_*` files — that is expected, not a failure. `amip_pi_forcing`
retains 12 coupler files (`A_SST_OpenIFS_*.nc`, `A_Ice_frac_OpenIFS_*.nc`); the rest have
been cleaned out. **None of them is a tuning lever and none should appear in an evaluation.**

*Practical note:* the evaluation scripts glob on `atm_remapped_1m_<var>_1m_<year>.nc`, so
these directories are skipped automatically with an "incomplete, skipped" warning if ever
added to a `RUNS` list. The warning is the intended behaviour, not an error to chase.

*Regenerate this inventory:*
```bash
cd /work/bb1469/a270092/runtime/oifsamip-cy48
for d in amip_*/; do
  printf "%-26s %s\n" "$d" "$(ls ${d}outdata/oifs/atm_remapped_1m_2t_1m_*.nc 2>/dev/null | wc -l)"
done
```

### ⚠ Round 13 results (2026-08-02): PREDICTION FALSIFIED — the snow lever is rejected

**The prediction was `+0.2` to `+0.7` K on Siberian JJA T2m. The measurement is `+0.020` K,
`t = 0.17`, pure noise.** Recorded as a clean failure. Worse, the lever inflicts the largest
winter cold bias in the campaign: **DJF −1.233 K**, against a ±0.588 K threshold — beating
B5's −0.720, which was rejected for exactly this.

| | DJF | MAM | JJA | SON |
|---|---:|---:|---:|---:|
| H1 `snowcr30` | **−1.233*** | −0.473* | **+0.020** | +0.428 |
| H2 = G1 + snow | −0.973* | −0.342 | +0.662* | +0.863* |
| *(G1, for reference)* | +0.070 | −0.296 | +0.521* | +0.280 |

**Why it failed — the mechanism was pointed at the wrong month.** The lever worked exactly as
designed on albedo; it just does not bite where the reasoning assumed
(`scripts/analysis/monthly_lever_check.py`, new). H1 all-sky surface albedo change ×100:

| Jan | Feb | Mar | Apr | May | **Jun** | Jul | Aug | **Sep** | **Oct** | **Nov** | Dec |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −0.51* | −0.21 | −0.42* | −0.84* | −1.43* | **−0.95*** | −0.10* | −0.32* | **−3.08*** | **−8.16*** | **−4.66*** | −1.60* |

The response is concentrated in **October (−8.2), November (−4.7), September (−3.1)** — the
**snow-onset** season — and is small in June (−0.95). The error in the design reasoning is
now obvious: `ZCVS = min(1, d_cm · RQSNCR)` only responds where snow is **shallow**, and
snow is shallow when it is *accumulating in autumn*, not when it is melting in June. By June
the Siberian snowpack is either still deep enough to saturate `min(1, ·)` under both settings,
or already gone. **Shallow snow is an autumn phenomenon, not a melt phenomenon.**

Consistent with that, June/July/August T2m move by −0.076 / −0.067 / +0.205, none significant,
while **October T2m is +1.135 K** — the darker autumn surface does warm, in the season nobody
asked about.

**Why it cools winter so hard**, given the albedo also *falls* in Dec–Mar: it is not a
radiative effect. February and March are dark (control albedo is meaningless at that sun
angle) yet cool by −1.65 and −1.19 K. Reducing the snow-covered fraction moves area from the
exposed/sheltered **snow tiles (5, 7)** onto the **vegetation tiles**, and the snow tile is
what buffers the surface against radiative cooling in polar night. The lever therefore trades
a small autumn warming for a large winter cooling through **tile fractions, not albedo**.

**Verdict: reject.** It fails the target it was built for and damages the coupled model's
original complaint. `RQSNCR` is reverted; the tree is back at as-released.

**What is still worth keeping from it:**
1. The June surface-albedo bias is **not** dominated by the snow-cover-fraction formulation.
   Whatever sets it — snow albedo decay, the vegetation masking of snow albedo, or the
   albedo of the snow-free surface itself — is still unidentified.
2. `RQSNCR` has enormous leverage on **autumn** albedo (−8 points in October). If an autumn
   or snow-onset bias is ever diagnosed, this is the knob.
3. H2 posts the campaign's best SO SW RMSE (**4.787**), best subpolar N Atl (**4.713**,
   −0.294) and best global T2m RMSE (**1.530**) — but all of that comes from D2b, and the
   winter damage disqualifies the combination. **G1 remains the configuration to carry
   forward.**

### ⭐⭐ Round 14 RESULTS (2026-08-03): G4 adopted — best of campaign, mechanism falsified

| run | setting | Siberia JJA | DJF | verdict |
|---|---|---:|---:|---|
| **G2** | `RVRSMIN(3,4)=500` | +0.336* | −0.325 | 64 % of G1's gain |
| **G3** | `RVRSMIN(3,4)=2000` | +0.876* | **−0.798*** | **rejected** — winter damage |
| **G4** | G1 + `RVRSMIN(9)=225` | **+0.952*** (t=7.68) | −0.385 | **ADOPTED** |

**G4 is the best configuration of the campaign** and closes **48 %** of the boreal bias
(G1 closed 26 %):

| metric | control | G1 | **G4** |
|---|---:|---:|---:|
| Siberia JJA T2m [°C] | 9.73 | 10.25 | **10.68** (+0.952) |
| Siberia sfc net SW | 153.78 | 159.31 | **161.88** (+8.10) |
| SO SW RMSE | 6.877 | 4.809 | **4.800** |
| subpolar N Atl SW RMSE | 5.007 | 4.872 | **4.738** |
| Nordic Seas SW RMSE | 9.058 | 9.358 ⚠ | **9.147** (+0.089) |
| tropics net TOA | 42.61 | 42.71 | 42.73 |
| global T2m RMSE | 1.579 | 1.553 | **1.543** |
| DJF / MAM / SON | — | +0.07/−0.30/+0.28 | −0.385 / −0.360 / +0.185, all within noise |

G4 also **largely repairs G1's Nordic Seas damage** (+0.300 → +0.089), which nothing
predicted and nothing yet explains — an unexplained improvement is as much a loose end as
an unexplained degradation.

**⚠ THE MECHANISM IS FALSIFIED.** The prediction on record was a measurable drop in
May–June snow water equivalent. Siberian land box SWE [mm], run × year ANOVA thresholds:

| | Mar (±3.53) | Apr (±3.85) | May (±5.25) | Jun (±2.72) |
|---|---:|---:|---:|---:|
| control | 121.98 | 131.07 | 97.12 | 15.39 |
| G1 | +0.88 | +0.37 | **+5.58*** | **+5.68*** |
| **G4** | −1.44 | −2.39 | **+2.75** | **+1.60** |

**G4's snow mass did not move.** T2m rose +0.952 K with May/June SWE unchanged within
noise, so the gain is the **plain sensible-heat route**, exactly as for F4 — *not* the
melt-timing route that motivated the run. The tundra lever stands on its independent
table-consistency justification and it works; the mechanism claimed for it does not.

**And F4 makes the snow bias worse.** G1 *significantly increases* May–June SWE (+5.58,
+5.68). The campaign's best boreal lever improves temperature while **aggravating** the
snow bias it was supposed to help. Leading suspect: cutting evapotranspiration also cuts
**sublimation**, a genuine snowpack sink. Untested.

**Predictions vs outcome:**

| prediction | outcome |
|---|---|
| G2 keeps >60 % of G1's gain if saturating | ✅ 64 % |
| G3 adds < +0.15 K beyond G1 if saturating | ❌ **+0.355** — not saturating |
| G4 gains +0.2…+0.5 K beyond G1 | ✅ +0.431 |
| G4 gain comes with a **May–June SWE drop** | ❌ **FALSIFIED** |

**`RVRSMIN` does not saturate — cap it at 1000.** 500/1000/2000 → +0.336/+0.521/+0.876;
increments per doubling +0.185 then **+0.355**, i.e. *accelerating*. There is no knee.
1000 is defensible only as **the largest value before the winter damage starts** — G3 at
2000 costs DJF −0.798, clearing ±0.588 and joining B5 and H2 in the "warms JJA, cools DJF"
list. This is the outcome that was flagged in advance as an argument *against* the
approach, and it stands: F4's magnitude is a fitted ramp, not a corrected bias.

**Two mechanisms for the June albedo bias are now dead**: snow *cover formulation*
(round 13) and snow *melt rate* (round 14). Untried candidates: snow-albedo decay
timescale; the canopy-masking tile-7 formulation (which the round-13 residual showed is
already doing most of the work); and the `RVVEGALB` table, never touched.

### ⭐⭐ Round 14 (2026-08-03): the June albedo bias is a SPRING MELT bias — and tundra

Round 13 established what the June bias is *not*. This is what it **is**, from output already
on disk plus ERA5 — no new integration needed. `albedo_decompose.py`, `albedo_decompose_prep.sh`.

**1. Against ERA5 (land-masked, June, `amip_presentday` vs 1990–2014):**

| June | model | ERA5 | Δ |
|---|---:|---:|---:|
| surface albedo `fal` | 0.2147 | 0.1730 | **+0.042** |
| snow cover fraction | 0.380 | 0.261 | **+0.119** |
| snow albedo `asn` | 0.658 | 0.759 | **−0.102** |

Both `f_snow` columns are the **same** HTESSEL formula applied to each dataset's snow mass, so
the gap is **snow AMOUNT, not the cover formula** — exactly why round 13 could not move it.
And the model's snow is **too dark, not too bright**. Brightening snow pushes the wrong way.

**2. Accumulation is right; the melt is a month late and 30 % too slow.** SWE [mm]:

| | model | ERA5 | Δ | δ(mod) | δ(ERA5) |
|---|---:|---:|---:|---:|---:|
| Feb | 106.5 | 97.6 | +9.0 | +16.7 | +15.3 |
| Mar | 122.1 | 114.5 | +7.7 | +15.6 | +16.9 |
| **Apr** | 131.6 | 112.8 | +18.7 | **+9.4** | **−1.6** |
| **May** | 98.1 | 65.3 | **+32.8** | **−33.5** | **−47.5** |
| Jun | 14.1 | 9.6 | +4.5 | −84.0 | −55.7 |
| Nov | 49.2 | 41.5 | +7.7 | +25.9 | +25.4 |

Accumulation increments track ERA5 almost exactly. **April: the model still gains while ERA5
has already peaked — the pack peaks a month late. May: melts 30 % too slowly.**

**3. May is albedo, not cloud — and it is a runaway.** vs CERES:

| | SWnet mod | SWnet CERES | Δ | SW↓ mod | SW↓ CERES | Δ |
|---|---:|---:|---:|---:|---:|---:|
| Apr | 87.8 | 86.6 | +1.2 | 176.9 | 187.7 | −10.8 |
| **May** | 122.5 | 136.4 | **−13.9** | 220.7 | 224.2 | **−3.5** |
| Jun | 167.6 | 189.7 | −22.2 | 215.7 | 228.2 | −12.5 |
| Jul | 165.9 | 180.4 | −14.5 | 193.8 | 208.7 | −14.9 |

May's downwelling is nearly right while the net is −13.9 → **almost pure surface albedo**.
June is half cloud / half albedo; July is essentially all cloud. Loop: late melt → more May
snow → higher albedo → −13.9 W/m² → less melt energy → snow survives into June. This is why
the model melts slower *despite* darker snow: extent beats brightness.

**4. ⭐ Tundra is the largest cover type in the box, and F4 misses it.** From the model's own
`tvh/tvl/cvh/cvl`, area-weighted over land:

| | type | cover | |
|---|---|---:|---|
| low | **9 Tundra** | **25.6 %** | **F4 does not touch this** |
| high | 4 Decid Needleleaf | 19.2 % | F4 |
| high | 18 Mixed Forest | 5.2 % | |
| low | 13 Bogs/Marshes | 3.9 % | |
| high | 3 Evergrn Needleleaf | 3.4 % | F4 |

**`RVRSMIN(9) = 80 s/m` is the lowest value of any vegetated type in HTESSEL** — below crops,
short grass and tall grass (all 100). The model has arctic tundra transpiring more freely than
tropical grassland. Its closest analogues, evergreen and deciduous shrubs, are both **225**;
co-occurring bogs/marshes are 240. This is a **table-consistency** argument, checkable against
the table and independent of our bias — unlike F4's unanchored 4× excursion.

**Round 14 runs — all namelist-only, one binary, no rebuild, launched together:**

| run | setting | prediction |
|---|---|---|
| **G2** | `RVRSMIN(3,4)=500` | if saturating, keeps >60 % of G1's gain (≥ +0.31 K) |
| **G3** | `RVRSMIN(3,4)=2000` | if saturating, adds < +0.15 K beyond G1 |
| **G4** | G1 + `RVRSMIN(9)=225` | **+0.2 to +0.5 K beyond G1, concentrated in JUNE**, with a measurable drop in May/June SWE |

**G4's falsifier: if T2m rises with no change in May–June snow mass, the melt mechanism is
unsupported** and the gain is just the plain sensible-heat route. G2/G3 falsifier: if the
response is linear in `RVRSMIN` we are riding a ramp with no natural stopping point, which is
an argument *against* the approach, not for a bigger number.

### Tuning without rebuilding: `NAMSURFTUNE` (2026-08-02)

The ECMWF `surf` library ships no namelist, so every HTESSEL lever was a source edit plus a
rebuild — serialising experiments against the one model tree, and putting AWI tuning inside
upstream files where it collides with EC-Earth at every merge.

`surf/module/surfece.F90` is already an AWI module *inside* surf with its own namelist, so the
overrides live there: `ECE_TUNE_RVRSMIN`, `_RVLAI`, `_RVCOV`, `_RVZ0H`, `_RVLAMSK` (0:20 by
vegetation type) and scalar `ECE_TUNE_RQSNCR` — the whole F-series, H-series and skin levers.
Applied by `SURFECE_APPLY_TUNING` from `susurf.F90` **after** `SUSURF_CTL` fills the tables;
nothing in setup derives from these entries, so a late override equals editing in place.
Sentinel defaults mean a run setting nothing is bit-identical to untouched code, and every
override is logged to `NULOUT`.

Own `&NAMSURFTUNE` group, **not** `NAMECECFG` — that one is read **twice** from the same
`fort.4`, by `ECE_CONFIG` in `arpifs/ecearth.F90` and by `SURFECE_CONFIG` in surf, and a
Fortran namelist read aborts with *invalid reference to variable* on any name the reading
module does not declare. Putting the tune entries in `NAMECECFG` killed the arpifs read at
`su0yoma.F90:152`. **Found by a 1-day test run** — the argument for always doing one.
Exactly **one line** changes in an ECMWF file. oifs-48r1 commit `1004cba`.

### Round 13 design (2026-08-02): the H-series — snow cover fraction, the June lever

*Written before the runs finished, per the rule learned in round 12.*

**Why snow.** The albedo investigation decomposed the 12.5 W/m² Siberian JJA surface-SW
deficit into ~7.0 W/m² cloud and ~5.5 W/m² **surface albedo**, and the albedo half is
concentrated almost entirely in **June** (albedo bias +0.086; July and August near-perfect).
That is a snow-*melt-timing* signature, not a snow-*brightness* one: the model holds a bright
surface into early summer and then agrees with observations once the snow is gone.

Before touching snow albedo the question "is it too cold because snow is too bright, or too
bright because too much snow because it is too cold?" was tested directly. **Winter cold is
not albedo-driven** — February–April is simultaneously cold *and* dark in the model, which
rules out the albedo→cold direction for those months. That leaves melt timing in the growing
season as the part worth fixing, which is also the only part that matters for the forest.

**The lever.** Snow cover fraction in `surfbc_ctl_mod.F90:317-322` is

    ZCVS = MIN(1, snow_depth_cm * RQSNCR)

with `RQSNCR = 1/10` under `LESN09=T` (`sussoil_mod.F90:157`). A grid box is therefore
declared **100 % snow-covered at just 10 cm of snow**, and linearly below that — with no
dependence on vegetation height or orographic roughness. For boreal forest this is
physically wrong in an obvious way: 10 cm of snow does not hide a spruce canopy, and the
observed masking depth over forest is several times larger. H-series sets `RQSNCR = 1/30`.

*Why this is safe for the ice sheets:* Antarctica and Greenland carry snow far deeper than
30 cm, so `MIN(1, ...)` saturates at 1 for any value in this range and they are untouched.
Only **shallow, marginal, melting** snow responds — which is exactly the target.

**Runs.** `H1` = `RQSNCR` alone against control (clean attribution); `H2` = `RQSNCR` on top
of G1's F4 source change and D2b namelist (the deployable stack). `H2 − G1` measures the
interaction directly against `H1 − control`.

**Falsifiable prediction, on record before the results:** H1 recovers a good part of June's
13.4 W/m² albedo term, worth ~+0.2 to +0.7 K on Siberian JJA T2m, and does **little in July
and August**, where the surface is already snow-free and the residual deficit is cloud.
*If it warms July as much as June, the mechanism is not the one claimed here* and the result
should be treated as an accidental global brightening, not a melt-timing fix.

> **⚠ FALSIFIED — see the results section above.** JJA came out at **+0.020 K** (t = 0.17).
> The failure mode was not the one anticipated: it did not warm July as much as June, it
> failed to warm *any* summer month, because the albedo response landed in **September–
> November** instead. Shallow snow is an autumn phenomenon, not a melt phenomenon.

Guardrails to check besides the usual set: DJF and MAM in `seasonal_by_run.py` (a snow change
is a cold-season change first), NH−SH albedo, and the Nordic Seas SW RMSE that G1 degraded.

### ⭐ G1 = F4 + D2b results (2026-08-02): the best configuration of the campaign

G1 combines the two levers whose mechanisms and geographies are disjoint — F4 (`RVRSMIN`
250→1000, boreal stomatal resistance) and D2b (`RCL_INPSEA=0.2`, `RCL_INPPMIN=700 hPa`,
ocean ice-nuclei scaling). It is the first pairing that improves **both** targets at once.

| metric | control | G1 | Δ | target |
|---|---:|---:|---:|---:|
| **SO SW RMSE** (priority) | 6.877 | **4.809** | −2.067 | — |
| SO TOA SW CRE [W/m²] | −60.29 | −63.13 | −2.84 | −68.14 |
| SO cloud area [%] | 83.07 | 83.59 | +0.52 | 89.72 |
| **Siberia JJA T2m [°C]** | 9.73 | **10.25** | **+0.521** (t=4.22) | ≈12.2 |
| Siberia sfc net SW [W/m²] | 153.78 | 159.31 | +5.54 | 166.26 |
| global net TOA [W/m²] | +0.64 | +0.45 | −0.195 | ~0 |
| tropics net TOA [W/m²] | 42.61 | 42.71 | +0.096 | 45.11 |
| global T2m RMSE [K] | 1.579 | 1.553 | −0.026 | — |
| Nordic Seas SW RMSE | 9.058 | 9.358 | **+0.300** | — |

**It superposes — measured, not predicted.** Every previous combination got superposition
wrong *in sign* (AB, ABB8). G1 is additive to within noise on every metric:

| | F4 | D2b | sum | G1 actual |
|---|---:|---:|---:|---:|
| Siberia JJA T2m | +0.749 | −0.222 | +0.527 | **+0.521** |
| SO SW RMSE | −0.217 | −1.914 | −2.131 | −2.067 |
| SO TOA SW CRE | −0.079 | −2.640 | −2.719 | −2.840 |
| global net TOA | +0.095 | −0.283 | −0.188 | −0.195 |

The likely reason additivity holds here and failed before: F4 acts on a land surface flux in
the boreal summer, D2b on ice nucleation over ocean below 700 hPa. AB and ABB8 combined levers
that both modified the same cloud scheme in the same regime.

**The one cost is the Nordic Seas**, +0.300 SW RMSE against a ±0.052 threshold — worse than
additive (F4 −0.074 + D2b +0.258 = +0.184). Nordic RMSE has never responded to anything else
in the campaign, so this is the first lever to move it, in the wrong direction. It is a small
box and not the deep-water priority, but it should not be allowed to grow further.

**Remaining gap:** G1 closes 26 % of the boreal bias (0.52 of ~2.0 K) and 36 % of the SO CRE
gap. The SO **cloud-area** deficit of ~6 pp is still untouched by every lever ever tried.

### Seasonal audit (2026-08-02): the noise floor is not the same in every season

`scripts/analysis/seasonal_by_run.py` (new) runs the same run×year ANOVA as `noise_floor.py`
but per season. This exists because the campaign evaluated JJA only for eleven rounds, and the
coupled model's original complaint is a **cold-season** bias.

Siberian T2m detection thresholds differ by a factor of 2.4 across the year:

| season | control | sd(eps) | 95 % threshold |
|---|---:|---:|---:|
| DJF | −29.35 °C | 1.407 K | **±0.588 K** |
| MAM | −9.28 °C | 0.923 K | ±0.386 K |
| JJA | +9.73 °C | 0.580 K | ±0.242 K |
| SON | −10.41 °C | 1.031 K | ±0.431 K |

**Applying the JJA threshold to a winter delta overstates significance by 2.4×, and I did
exactly that** in a first pass, briefly recording D2b (−0.455) and B3 (−0.486) as
winter-damaged levers. Both are **within** the DJF noise floor. Corrected findings:

- **B5 `capdcycl0` is the only lever with genuine winter damage** (DJF −0.720, clears ±0.588),
  and the only one that warms JJA significantly while cooling DJF significantly. Its rejection
  stands, now on a properly-thresholded basis.
- **D2b's real seasonal cost is spring, not winter**: MAM −0.467 clears ±0.386; its DJF −0.455
  does not clear ±0.588.
- **F4 is seasonally clean** — DJF +0.061, MAM +0.047, SON +0.110, all within noise, with the
  entire signal in JJA (+0.749). This is what a well-targeted lever looks like and is a further
  argument for F4 over the cloud levers.
- **G1 is seasonally clean**: DJF +0.070, MAM −0.296, SON +0.280, none significant. The MAM
  cooling inherited from D2b is diluted below its threshold.
- Also newly significant on this test: A1a cools MAM/JJA and warms DJF/SON by large amounts —
  it is not a boreal lever, it is a global cloud change. B6, B7 and F3 all cool MAM
  significantly; F5 inherits F3's MAM −0.507.

*Run it with every round from now on.* Reporting JJA alone is how B5 survived eleven rounds.

### ⭐ Surface albedo (2026-07-31): half the Siberian SW deficit is not cloud at all

Every boreal lever so far has attacked cloud. Decomposing the 12.5 W/m² Siberian JJA surface
SW deficit against CERES on the identical box and mask shows cloud is only half of it
(`scripts/analysis/albedo_by_region.py`, `albedo_vs_t2m_bias.py`):

| term | W/m² | mechanism |
|---|---:|---|
| too little SW reaching the surface | **7.0** | excess cloud — what F4 attacks |
| **surface albedo too high** | **5.5** | **independent of cloud** |

Model Siberian JJA albedo **0.1753** vs CERES **0.1461**. F4 does not touch it (0.1749): it
raises SW *down* by thinning cloud while the albedo error passes straight through.

#### It is not a global albedo error — temperate and tropical land is right

| region | model | CERES | diff | W/m² |
|---|---:|---:|---:|---:|
| **Siberian tundra** | 0.3364 | 0.2427 | **+0.094** | **−16.0** |
| Canada boreal | 0.1697 | 0.1382 | +0.032 | −6.1 |
| Siberia boreal | 0.1753 | 0.1461 | +0.029 | −5.5 |
| Sahara | 0.3811 | 0.3616 | +0.020 | −4.8 |
| **Fennoscandia** | 0.1157 | 0.1141 | **+0.002** | −0.3 |
| Europe / Amazon / Plains / steppe / Australia | — | — | **±0.004** | ~0 |

Temperate and tropical vegetation is essentially perfect, so the albedo *scheme* is not
broken and a global fix would be wrong. **Sahara +0.020 is a separate, real finding** worth
its own look. Note Fennoscandia — same needleleaf type as Siberia — is correct, which already
argues against a needleleaf-albedo error.

#### Latitude split: it is not the forest

| band | model | CERES | diff |
|---|---:|---:|---:|
| 55–60N | 0.1369 | 0.1257 | +0.011 |
| 60–65N | 0.1377 | 0.1272 | +0.011 |
| 65–70N | 0.1738 | 0.1429 | +0.031 |
| 70–75N | 0.3026 | 0.2152 | **+0.087** |

#### How the albedo is actually built — and why two hypotheses died

All configurations (AMIP, coupled round-09, **and the CMIP7 production piControl**) run
`LRDALB = .false.` with `NALBEDOSCHEME = 3`. That **discards the MODIS albedo sitting in the
input file** (params 15–18, 174) and constructs albedo every step (`surfbc_ctl_mod.F90:512`):

```
alpha = RVVEGALB(tvl)*cvl + RVVEGALB(tvh)*cvh + alpha_soil*(1-cvl-cvh)
PALBF = 0.45976*PALUVD + 0.54024*PALNID
```

*Hypothesis 1 — the soil-albedo field dominates tundra. **WRONG.*** The bare fraction is only
**0.01–0.03** in every Siberian band, so the soil term is negligible.

*Hypothesis 2 — the soil fields are corrupt. **Cosmetic only.*** GRIB 117–120 declare `N=128`
with a full 256-row `pl` array while carrying 40 320 values (the O96 count) — a genuinely
inconsistent header that **breaks `cdo` and `grib_get_data`** on those fields. But the values
are O96-ordered and correct: soil albedo is zero over 99.2 % of ocean cells and non-zero over
98 % of land, `corr(soil, lsm) = +0.894`. **Report upstream as a metadata bug; it is not our
bias.**

The blend at 70–75N gives only **0.158**, while the model runs at **0.303**. Vegetation and
soil cannot produce that — the extra ~0.145 is **snow**.

#### The real signature: snow-melt timing, not albedo magnitude

| band | Apr | May | **Jun** | Jul | Aug | **Sep** |
|---|---:|---:|---:|---:|---:|---:|
| 55–65N | −0.016 | +0.099 | +0.035 | −0.006 | +0.003 | +0.042 |
| 65–70N | **−0.061** | +0.029 | **+0.084** | +0.004 | +0.005 | **+0.082** |
| 70–75N | **−0.067** | −0.043 | **+0.116** | +0.092 | +0.055 | **+0.101** |

The model is **too DARK in Apr** (−0.06 to −0.07), **near-perfect in Jul/Aug** when the surface
is snow-free (+0.004), and **too bright in Jun and Sep**. So `RVVEGALB` is essentially right;
the snow season is **too long at both ends**.

#### Cause and effect: cold is upstream, but the June snow is a real defect

| | Feb | Mar | Apr | May | **Jun** | Jul |
|---|---:|---:|---:|---:|---:|---:|
| T2m bias | −1.70 | −1.41 | −0.08 | −0.34 | −1.35 | −2.35 |
| albedo bias | **−0.084** | **−0.089** | −0.063 | +0.005 | **+0.086** | +0.030 |

**Feb–Apr the model is too cold *and* too dark simultaneously.** If bright snow caused the
cold, the signs would agree. They do not — so the winter cold is **not** albedo-driven and
must come from elsewhere (the global tropospheric cold bias, longwave, or advection). The
albedo bias only turns positive in June, *after* the temperature crosses freezing.

But the June snow is not merely inherited: **June is +6.2 °C with 213 W/m² of insolation and
0.020 m w.e. still lying.** Snow should not survive that. Implied snow-cover fraction is
**~17 %, where CERES implies ~2 %** — a very thin layer credited with far too much area, which
points at the **snow-cover-fraction formulation**, not snow albedo (too *low* in April) and not
melt energy (ample).

#### Growing season: June is albedo, July is cloud

| Siberia 55–75N | Jun | Jul | Aug | JJA |
|---|---:|---:|---:|---:|
| cloud | 6.9 | **9.9** | 4.4 | 7.0 |
| **albedo/snow** | **13.4** | 1.7 | 1.3 | 5.5 |
| total | 20.2 | 11.6 | 5.7 | 12.5 |

Almost cleanly separated in time. The residual snow/albedo term is **5.5 W/m² on the JJA
mean, essentially all in June** — worth **~0.2–0.7 K** depending on which sensitivity applies
(0.033 K/W/m² from the spatial regression, 0.136 from F4's own response). Potentially another
F4-sized gain, and **orthogonal to F4**, which fixes SW-down while this fixes the absorbed
fraction.

#### But albedo does not explain the land cold bias

Insolation-weighted over land 60S–75N (excl. Greenland): mean albedo-induced SW loss
**−3.00 W/m²**, mean T2m bias **−1.35 K**, spatial correlation **r = +0.145**, slope
**+0.033 K per W/m²**. Decisively, **cells where the model is too *dark* are still −1.16 K
cold**. So a ~1.2 K land-wide cold bias has **no albedo signature at all**, consistent with the
global tropospheric cold bias of 0.7–2.9 K found in the vertical profiles. Fixing albedo is
worth doing; it is not the root cause.

*Caveats:* CERES surface albedo over snow at high latitude carries real retrieval uncertainty,
and the negative Apr/May bias is odd enough that ERA5 `fal` should be used as a cross-check
before trusting the amplitude. In a coupled run with LPJG this interacts with the
vegetation–albedo feedback that AMIP deliberately switches off.

### ⭐ Round 12 results (2026-07-31): F4 is the best boreal lever of the campaign

Threshold ±0.242 K at 44 yr. **F4 nearly doubles the previous best and costs no guardrail.**

| run | ΔSiberia JJA T2m | t | verdict |
|---|---:|---:|---|
| **F4** `RVRSMIN` 250→1000 | **+0.749** | +6.06 | **significant — best ever** |
| **F5** all four | **+0.746** | +6.03 | significant, ≡ F4 |
| B5 (previous best) | +0.407 | +3.29 | significant |
| F2 `RVLAI` 5→3 | +0.321 | +2.59 | significant |
| F1 `RVZ0H` z0m/10 | +0.185 | +1.49 | marginal |
| F3 `RVCOV` 0.9→0.7 | ≈0 | — | noise |

| | control | B5 | **F4** | F5 | target |
|---|---:|---:|---:|---:|---:|
| Siberia JJA T2m | 9.73 | 10.13 | **10.48** | 10.47 | ~12.2 |
| Siberia sfc SW | 153.78 | 156.26 | **159.33** | 160.04 | 166.26 |
| Siberia cloud % | 78.14 | 77.02 | **75.89** | 75.52 | 69.59 |
| global net TOA | 0.64 | **−0.36** | 0.74 | 0.69 | ~0 |
| **tropics net TOA** | 42.61 | **40.67** | **42.68** | 42.67 | 45.11 |
| T2m RMSE vs ERA5 | 1.58 | 1.56 | **1.54** | 1.56 | |
| Bowen ratio | 0.44 | 0.42 | **0.60** | 0.57 | 0.5–1.5 obs |

**F4 costs the tropics nothing** (42.68 vs control 42.61) where B5 wrecked them (40.67). That
is what vegetation-type indexing bought, and it worked as designed. Bowen moves 0.44 → 0.60,
into the observed range **without overshooting** — the risk flagged in advance did not
materialise.

**F5 ≡ F4 to 0.003 K.** Against a naive sum of ~1.25 K, adding F1+F2+F3 on top of F4
contributes *nothing*. The advance prediction ("F5 ≪ the sum, closer to the largest single
lever") held, and more strongly than expected: near-total saturation of one shared
latent-heat pathway. **Use F4 alone.**

**Prediction that failed:** F1 (`RVZ0H`) was called strongest a priori because kB⁻¹ = 0 is the
largest departure from observation. It is marginal (+0.185) and **barely moved Bowen**
(0.44→0.43). "Largest departure from observation" is not a reliable guide to leverage.

**B5 carries a hidden −0.72 K DJF cooling** (F4: +0.06). Invisible for the whole campaign
because only JJA was ever evaluated. Given the coupled model's original complaint is a
*cold-season* bias, that is disqualifying for B5. **DJF now belongs in the standard
guardrails.**

**Is F4 fixing a cause or a symptom?** Two tests, both partly exonerating:
* *Seasonally selective* — DJF +0.061, MAM +0.047, **JJA +0.749**, SON +0.110. Exactly the
  signature a transpiration mechanism must have (stomata shut, canopy snow-covered in winter),
  so it acts through the mechanism claimed.
* *Vertically targeted* — warms 1000/925/850 hPa by +0.83/+0.84/+0.83, closing **56–83 % of
  the Siberia-specific** layer bias, and only +0.53 aloft. It fixes the layer that is
  specifically Siberian rather than smearing warmth through the column.

What remains true: `RVRSMIN` 250→1000 is an **unanchored 4× excursion** from what ECMWF
ships, and the vertical profiles show it removes moisture that is already deficient.

### Round 11 results (2026-07-30): D2b adopted for the Southern Ocean, D1 falsified

| | SO SW CRE | SO SW RMSE | Sib JJA T2m | global TOA | tropics |
|---|---:|---:|---:|---:|---:|
| control | −60.29 | 6.88 | 9.73 | 0.64 | 42.61 |
| **piCTRL 1850** | −60.65 | 6.83 | 9.74 | **0.79** | 42.69 |
| D1 `RCAPDCYCL`=4 | −60.38 | 6.51 | **9.56** | 0.61 | 42.53 |
| D2a INP, no gate | **−65.13** | 5.07 | 9.56 | −0.28 | 42.23 |
| **D2b INP + p700** | −62.93 | **4.96** | 9.51 | 0.36 | **42.70** |

**D2b adopted.** Best SO SW RMSE of any run (6.88 → 4.96, −28 %) with **zero tropical cost**
(+0.09), beating A1b (5.54) and even D2a (5.07) despite a smaller mean-CRE change — it
improves the spatial *pattern*, which is what the priority metric is for.

**D1 falsified cleanly.** Predicted to recover part of B5's +0.407 K by reformulating the land
CAPE closure; delivered **−0.170 K**, the wrong sign. So B5's gain requires *removing* the
land correction, not rescaling it — the effect is not about which quantity the closure scales
on. The advance prediction made this unambiguous.

**piCTRL: the energy target is +0.79, not +0.64** — slightly *larger*, so A1b does not
overshoot. And **boreal T2m differs by only +0.008 K**, so the dead-namelist forcing bug never
affected the boreal target and every round-10 boreal conclusion stands.

**G1 = F4 + D2b is running** — the first pairing where the two targets do not conflict, since
they are disjoint in both process and geography (stomatal conductance on vegetation types 3/4
versus ice nucleation over ocean below 700 hPa). Naive sum: boreal +0.53, SO RMSE ~4.74,
tropics ~42.8. **A prediction, not a measurement** — AB and ABB8 both got superposition wrong
in sign, and F5 saturated almost completely.

*Infrastructure:* `eval_round10_A.py` now caches per run (JSON in `.eval_cache/`) and
parallelises the cold path — **2m09s → 7.7s, output verified byte-identical**. Adding one run
costs one run's work. `table.py` builds tables from the cache; the old fixed-width text parser
silently mis-aligned columns when new run labels had different widths.

### Round 12 design (2026-07-31): the F-series — boreal surface exchange

Round 11 left the boreal unsolved: B5 (+0.407 K) is still the only significant lever and it
costs **−1.94 W/m² in the tropics**, which are already 2.5 below CERES. D1 was built to
separate B5's boreal gain from that cost and **failed** (−0.170 K, wrong sign). So round 12
changes target: instead of convection or cloud, it attacks **boreal surface exchange**.

**One observational anchor for all five runs.** The measured Siberian JJA Bowen ratio is
**0.43**, against ~0.5–1.5 observed for boreal conifer — the model puts too much of the
surface energy into latent heat. Each lever raises Bowen by an independent route, so this is
one hypothesis tested four ways rather than four unrelated pokes.

**All are indexed on vegetation types 3 and 4** (evergreen/deciduous needleleaf), so types
5/6/18 (broadleaf, mixed) keep their defaults and the tropics and temperate forests are
untouched **by construction**. That is precisely the property B5 lacks, and the reason an
F-lever could beat it even at similar boreal magnitude.

| run | change (`susveg_mod.F90`) | basis | binary |
|---|---|---|---|
| **F1** | `RVZ0H(3,4)` = `RVZ0M/10` | see below — largest departure from observation | `feb42a8c3985` |
| **F2** | `RVLAI(3,4)` 5.0 → 3.0 | observed boreal conifer LAI is 2–4, not 5 | `12ed4f92e7f5` |
| **F3** | `RVCOV(3,4)` 0.9 → 0.7 | taiga is often open woodland | `03d49cb864d0` |
| **F4** | `RVRSMIN(3,4)` 250 → 1000 | expA's 250→500 gave +0.232 K; tests saturation | `fafde82f9d91` |
| **F5** | all four together | **measures** the superposition | `2b6aaac9da20` |

None had been tried before: `RVZ0H`, `RVLAI` and `RVCOV` had **zero** prior mentions in this
logbook, and `RVRSMIN` only at 500 (expA). Five distinct md5s, each verified against what its
experiment staged, with the revert verified between every build so each of F1–F4 isolates one
parameter and F5 is exactly their sum.

#### F1: forests use `z0h = z0m`, everything else uses `z0m/100`

The find that motivated the round. In `susveg_mod.F90`, **every closed-forest type (3, 4, 5,
6, 18) sets the roughness length for heat equal to that for momentum**, while every other
land type uses `z0m/100` (water/ice use `z0m/10`). Even "Interrupted Forest" (19) gets /100.
So it is systematic and deliberate, not a typo — but observations give
kB⁻¹ = ln(z0m/z0h) ≈ 2 for forests (z0h ≈ z0m/7), i.e. the model uses **kB⁻¹ = 0** and
parameterises forest heat *and moisture* exchange as far more efficient than observed.
Raising the resistance forces the surface to warm to shed the same energy and cuts
evaporation — both raise Bowen.

*Honest caveat on sign:* for **skin** temperature the sign is unambiguous (warmer). For **2 m**
temperature it is not — reducing `z0h` steepens the near-surface gradient and T2m is diagnosed
within it, so part of the skin warming may not reach 2 m. The moisture side reinforces
warming. Net warming is the **prediction**, not a certainty.

#### F5: measure the superposition, do not infer it

Predicting superposition has failed **twice, in sign**: AB (A1b + B2) predicted +0.21 K and
measured −0.053; ABB8's parts summed to +0.76 K and measured −0.167. F5 costs one build and
one job and runs in parallel, so the superposition arrives *with* the individual sensitivities
rather than nine hours later.

**Prediction on record: F5 ≪ F1+F2+F3+F4**, probably closer to the largest single lever than
to their sum. Unlike AB — which combined genuinely different processes — F1–F4 all act through
the **same pathway**, suppressing latent heat flux, and latent heat has a floor at zero, so
saturation is close to guaranteed. If F5 comes out near-additive, the four are working through
partly separate channels (e.g. F1 altering boundary-layer structure rather than only
evaporation), which is itself worth knowing.

**Risk named in advance:** all four at full strength may drive Bowen *past* the observed
0.5–1.5 range into a dry bias. A large F5 is therefore **not automatically good** — check
Bowen, soil moisture and T2m RMSE, not just the Siberian JJA mean. With F1–F4 individually
plus one superposition point we can interpolate down to a sensible amplitude instead of
guessing.

#### What the vertical profiles say about this round

Measured *after* the runs were launched (see the section above), and it cuts both ways:

* **Supports the framing.** The Siberia-specific cold bias is confined below 850 hPa
  (−1.0 to −1.5 K beyond the global bias, ~0 at 700–300 hPa). That is exactly the signature a
  surface-exchange problem should have, so the F-levers are aimed at the right layer.
* **Cautions the mechanism.** Siberian RH is +5.6 to +8.2 % too high while `q` is **too low**.
  The excess cloud is *thermally* driven — cold air saturating at lower moisture — not
  moisture-driven. The F-series removes moisture that is already deficient, so **a large
  F-response must not be read as confirmation that we found the cause.**

Also relevant: above 700 hPa Siberia merely shares a **global tropospheric cold bias of
0.7–2.9 K**. No boreal lever will touch that, and it is a separate problem the campaign has
not yet addressed.

### Vertical structure vs ERA5 (2026-07-31): where the biases actually live

Everything measured until now was at the surface or TOA, which says a bias exists but not
where it originates. `scripts/analysis/vertical_profiles_prep.sh` + `vertical_profiles.py`
compare model and ERA5 on the model's 19 pressure levels.

*Setup:* ERA5 monthly pressure levels come from the **DKRZ pool**
(`/pool/data/ERA5/E5/pl/an/1M/`, params 130 T, 133 q, 157 RH) — **no download needed**. All
19 model levels exist *exactly* in ERA5's 37, so `-sellevel` matches with no interpolation.
ERA5 is regridded onto the model grid so the *same* land mask applies to both. The model
writes **monthly** `pl` output (36 MB/yr/var) as well as 6-hourly (4.6 GB) — use the monthly.
Primary comparison is `amip_presentday` vs ERA5 1990–2014, which is **period-clean**.

#### ⚠ The model's pressure-level RH field is BROKEN

`r` on pressure levels is **identically zero in every run checked** — 0 nonzero of 17,510,400
values in both `amip_presentday` (1995) and `amip_pi_base` (1900) — while `q` on the same
levels is fine. The XIOS `r` output is dead and **must not be used**; it silently yields
garbage rather than failing. RH below is computed from `t` and `q` via Bolton's formula
applied identically to both datasets, so the saturation convention cancels in the difference.

#### The Siberia-specific cold bias is confined below 850 hPa

| hPa | Siberia | Global | **Siberia-specific** |
|---:|---:|---:|---:|
| 1000 | −1.90 | −0.73 | **−1.17** |
| 925 | −2.17 | −0.65 | **−1.52** |
| 850 | −1.99 | −0.99 | **−1.00** |
| 700 | −1.22 | −1.15 | −0.07 |
| 500 | −1.67 | −1.49 | −0.18 |
| 300 | −2.32 | −2.22 | −0.10 |
| 200 | −4.92 | −2.87 | −2.05 |

**This separates two problems that had been conflated.** Above 700 hPa Siberia merely shares
a **global tropospheric cold bias of 0.7–2.9 K** (peaking near the tropopause) that no boreal
lever will touch. The boreal-specific excess lives entirely in the bottom three levels —
exactly the signature a surface-exchange problem should have, which **validates the F-series
framing** (`RVZ0H`/`RVLAI`/`RVCOV`/`RVRSMIN` all act on surface exchange).

#### The excess boreal cloud is thermally driven, not moisture-driven

Siberian JJA RH is **+5.6 to +8.2 % too high**, but `q` is **too LOW** (−0.2 to −0.3 g/kg).
The air is not moist; it is cold, so it saturates at lower moisture. A 2 K cold bias cuts
`q_sat` by ~14 % (Clausius–Clapeyron) while `q` is only ~6 % low, which accounts for the RH
excess almost exactly.

**Caution for the F-series.** Those levers work by *removing moisture*, which will reduce
cloud — but moisture is already deficient and the loop reads cold → cloud → less SW → colder.
Cutting an already-deficient moisture supply may be treating a symptom. A large F-response
should therefore **not** be read as confirmation that we found the cause.

#### Independent support for the D2 (INP) hypothesis

The Southern Ocean has RH **+5.8 to +6.6 % too high through the mixed-phase layer**
(850–500 hPa) while reflecting **too little** SW (CRE −60.3 vs CERES −68.1). The humidity is
present; the liquid is not. That is precisely the supercooled-liquid deficit the
ice-nucleation hypothesis predicts, seen from a completely independent direction.

*Qualification:* the SO RH excess actually **peaks at 300 hPa (+10 %)**, above the layer
D2b's 700 hPa gate targets. So the gate's success may owe more to sparing the tropics than
to selectively hitting SO low cloud, and the low-cloud framing in the round-11 notes is
weaker than the D2a-vs-D2b result alone suggested.

*Not comparable:* the model writes no cloud variables on pressure levels (`cc`/`clwc`/`ciwc`
absent), so vertical cloud structure cannot be compared directly — RH is the only proxy.
Adding those three to the `file_def` would be a small change worth making for future runs.

### Round 11 design (2026-07-30): two physically-motivated levers

Round 10 ended with two structural walls (see the combination arithmetic above): the SO is
reachable only through mixed-phase overlap that cools Siberia, and Siberia's levers sum to
+1.06 K against a +1.3–2.0 K need without adding. Round 11 attacks the *mechanisms* instead
of trading parameters.

**Two proposals were withdrawn first, for being tuning-to-fit rather than physics.**

1. ~~Land/sea split on `RCL_OVERLAPLIQICE`.~~ **Withdrawn — anti-physical.** That parameter is
   the sub-grid **overlap fraction of supercooled liquid with ice** (range [0,1], default
   0.65, EC-Earth4 uses 0.1), entering `ZDEPOS = MAX(ZOVERLAP_LIQICE·ZA·(ZINEW−ZICE0),0)`.
   The code's own comment gives the intended selectivity: *"Reduce in shallow convection
   because assume SLW in active updraught is less overlapped with ice in less active part"* —
   i.e. **convective segregation lowers overlap**. But ocean mixed-phase cloud is
   predominantly *stratiform* (well mixed → HIGH overlap) and land cloud more convective
   (segregated → LOW overlap), which is the **opposite** of what the tuning wants. No
   physical story; it only helps us.
2. ~~`RCAPDCYCL` mode 3 = disable the CAPE diurnal-cycle correction over land.~~
   **Withdrawn — a regression.** That correction is Bechtold et al. (2014), which exists to
   fix the documented error of land convection triggering too early in the day. Switching it
   off is not a fix.

#### D1 — `RCAPDCYCL` mode 4: reformulate the land closure, do not delete it

**B5 was misdiagnosed.** `RCAPDCYCL` is a **mode selector**, not a magnitude
(`cumastrn.F90:773,777`):

| value | behaviour |
|---|---|
| 1.0 | correction over **land only**, scaled on `ZKHVFL` (surface kinematic virtual heat flux, `:464`) |
| 2.0 | **everywhere** — land via `ZCAPPBL·ZTAU/ZTAURES`, ocean via a wind-based `ZTAUPBL` |
| 0.0 | **off entirely** ← what B5 did |

So B5's +0.407 K did not tune anything down: it **deleted the correction globally**. Its
boreal gain came from the land branch, its **−1.94 W/m² tropical cost** from removing the
ocean branch. The two are separable.

**The physical argument** is in what the land formulations scale on. Mode 1 is proportional
to the *actual surface buoyancy flux* and **vanishes when that forcing is weak**; mode 2's
land branch uses boundary-layer CAPE instead. In boreal summer — weak insolation, shallow
BL, small surface fluxes — the buoyancy-flux version self-limits while the BL-CAPE version
need not. That is a real reason to think **mode 2's land branch is mis-scaled for
weak-forcing high-latitude convection**, and it predicts the sign we measured.

**D1 adds mode 4 = mode-1 land physics + mode-2 ocean branch verbatim.** Both branches are
existing code; no new physics is invented. Falsifiable prediction: recovers part of B5's
+0.407 K with little of its tropical cost.

#### D2 — ice nuclei, not overlap: `ZICENUCLEI` continental scaling

The SO supercooled-liquid deficit is an **ice-nucleation** problem, and the scheme has no
aerosol dependence at all:

```
ZICENUCLEI = 1000·EXP(12.96·ZSUPERSATICE − 0.639)      (cloudsc.F90:2516, :2618)
```

That is **Meyers et al. (1992)** — a function of ice supersaturation *only*, documented to
overestimate INP in clean environments, which is exactly the remote Southern Ocean.
DeMott et al. (2010/2015) make INP a function of aerosol number instead. Scaling
`ZICENUCLEI` down over ocean therefore has a **real mechanism** — INP sources are
continental (mineral dust, biological) — with the correct sign in both regions at once:
low INP over the SO → liquid survives → brighter; continental INP unchanged → boreal
untouched. Leverage is strong since `ZCVDS ∝ ZICENUCLEI^0.666`, so a factor 0.1 cuts
deposition to ~0.22.

**Why the "proper" version is not available.** `PNICE` (ice number concentration) is already
a `cloudsc` argument (`:265`) and already used for autoconversion (`:2855`), but it is filled
only when `NAERCLD > 0` (`callpar.F90:1155`, else zeroed at `:1166`); `NAERCLD` defaults to 0
(`sucldp.F90:597`). **`NAERCLD > 0` is ruled out on cost** — that configuration is several
times slower and would not finish CMIP7 in time.

**Why the dust climatology cannot be used either.** `onlinedust_v4_rmp.nc` is staged and
*opened* unconditionally, but its arrays (`POTSRC`, `SOILTYPE`, `CULT`) are allocated only by
`TM5M7_SRC_DUST_INIT` via `tm5m7_init.F90` — i.e. **only under m7**. `suecrad`'s `YDUSTCLIM`
is a different object (dust optical properties for the radiation aerosol climatology, not the
source map). So the file's contents are never populated in our configuration.

**Why MACv2-SP is not a usable INP proxy.** Its plumes are anthropogenic sulphate/BC/OC —
poor INP, and sulphate coating *deactivates* mineral dust INP. Worse, in a piControl
anthropogenic aerosol is near zero, so it would give almost no land/ocean contrast.

**Two physical caveats that shape the test — added after review.**

*The tropics will pay, and "ocean" is the wrong discriminator.* The scaling acts only inside
the mixed-phase window (`RTHOMO ≤ T < RTT−5`, i.e. 235.15–268.15 K, liquid present), which in
the tropics sits at ~500–300 hPa and **is** populated: deep convective updraughts carry
supercooled liquid well above the freezing level, and anvil detrainment feeds it. Reducing
INP there retains more liquid → brighter tropical mid-level cloud → tropical TOA pushed
**further below** CERES, where it is already 2.5 too low. Expect a cost of order A1a's
−1.07 W/m² or worse. And physically, tropical/subtropical oceans are **not** dust-remote
(Saharan dust over the tropical Atlantic, Asian dust over the N Pacific) — the SO is uniquely
far from sources. The distinguishing property is **remoteness**, which a local land mask
cannot encode.

*Sea salt sets an INP floor.* Sea salt itself is a poor INP (highly soluble, inefficient in
the immersion mode at these temperatures), but sea spray carries **marine biogenic organics**
— sea-surface microlayer material, phytoplankton exudates — which are a modest INP source
(Wilson et al. 2015; DeMott et al. 2016). So **the ocean factor must not approach zero:
~0.1–0.3 is defensible, not 0.01.** Sea salt's larger role is as CCN, which is moot here since
`PCCN` is zeroed at `NAERCLD=0` and droplet number is the constant `RCL_KK_CLOUD_NUM_SEA`.
Tempering caveat: marine biogenic INP peaks with the austral-summer bloom (DJF), exactly when
SO insolation and our SW bias peak — so the SO is *least* INP-starved in the season we care
most about, which argues against expecting A1a-sized leverage.

**D2 phase 1 (this round):** scale `ZICENUCLEI` by a continuous function of `PLSM` (land
fraction, already passed and already used land-selectively at `:2054/:2062/:2113/:2131`), with
a namelist factor so it can be scanned without rebuilding.

**D2 phase 2 (if phase 1 has leverage):** replace the land-fraction proxy with a
**tile-based two-term INP field** computed in `callpar` from `SURFL%ZFRTI` (already used at
`callpar.F90:2006`) and passed into `cloudsc`:

```
INP ∝ a·[bare soil (tile 8) × wind above threshold]   <- mineral dust
    + b·[low veg (tile 4) + high veg (tile 6)]        <- biological
```

Tile convention from `srfi_mod.F90:52`. This is better physics *and* gets the signs right
where a dust-only proxy would not: boreal forest is a weak dust source but a **strong
biological INP source in summer** (pollen, bacteria, fungal spores, litter), which is exactly
our season — while snow tiles (5/7) collapse the biological term in winter, leaving the
cold-season bias untouched. A dust-only proxy would *lower* boreal INP and brighten boreal
cloud, i.e. the wrong sign.

**Build note:** both D1 and D2 are source edits on the single `oifsamip-cy48` tree, so they
**serialise** — edit → `esm_master recomp-oifsamip-cy48` → verify binary md5 changed → submit
→ wait for staging → next. (`recomp-awiesm3-develop-is` rebuilds the *coupled* tree, which
our AMIP runs do not use; propagating these edits there is a separate step for the coupled
rounds.)

### 44-year results (2026-07-30): the boreal ranking is rebuilt

Detection threshold fell from ±0.89 K to **±0.240 K** (`sd(eps)`=0.575, dof=774).

| run | 4-yr Δ | **44-yr Δ** | t | verdict |
|---|---:|---:|---:|---|
| **B5** `RCAPDCYCL` 2→0 | +0.001 | **+0.407** | +3.32 | **SIGNIFICANT — best** |
| **B2** `RCLDIFF_CONVI` 25 | +0.290 | **+0.344** | +2.80 | **SIGNIFICANT** |
| **B3** `RCLDIFF` 1.5E-5 | +0.366 | **+0.306** | +2.50 | **SIGNIFICANT** |
| A1a ovl=0.10 | −1.149 | −0.749 | −6.11 | SIGNIFICANT (wrong way) |
| expA, C1 | +0.193, −0.213 | +0.232, +0.223 | ~1.9 | marginal |
| **B8** `RVLAMSK`=5 | **+0.502** | **−0.038** | −0.31 | **noise** |
| AB, B7, C2, E1, A2, A1c, B6 | — | ≤ +0.12 | — | noise |

**B8, the "best lever of the campaign", is noise.** The 07-29 retraction was correct.
**B5, dismissed as "inert" at +0.001, is the strongest lever.** The caution written then —
"unresolved, not proven inert" — is what preserved it.

**Guardrails at 44 yr — and a real conflict:**

| metric | control | B5 | B2 | B3 | A1b | AB | ABB8 |
|---|---:|---:|---:|---:|---:|---:|---:|
| SO SW CRE | −60.29 | −60.52 | −60.18 | −59.62 | **−62.46** | −62.34 | −62.17 |
| Sib JJA T2m | 9.73 | **10.13** | 10.07 | 10.03 | 9.59 | 9.68 | 9.71 |
| global net TOA | 0.64 | **−0.36** | 1.04 | 1.98 | **0.07** | 0.39 | 0.43 |
| tropics net TOA | 42.61 | **40.67** | 43.26 | 44.62 | 42.32 | 42.83 | 42.92 |

* **B5** wins boreal *and* improves |TOA| (0.64 → −0.36), but costs **−1.94 W/m² in the
  tropics**, already 2.5 too low. That is its price.
* **B2/B3 buy boreal with energy** (+0.39, +1.34 TOA). B3 is unusable.
* **A1b nearly perfects energy** (0.64→0.07) and owns the SO — but **cancels B2's boreal
  gain**: additive prediction +0.21, measured **−0.053**. The anti-synergy is now solid.
* **A1b and B5 both push TOA negative** (−0.58, −1.00), so stacking them would overshoot to
  ≈−0.94. The best boreal lever and the best energy lever are not freely combinable.

**Deep-water SW RMSE at 44 yr** (per-year RMSE, ANOVA): SO threshold ±0.019 — A1b −0.12,
AB −0.12, ABB8 −0.11 all significant. Subpolar N Atl ±0.024 — **B2 +0.029 significant
damage; AB −0.023 noise; ABB8 −0.026 a significant *improvement***, so the cancellation is
confirmed and slightly better than neutral. **Nordic Seas: nothing significant** (±0.053) —
A1b's −0.989 was never established. Global SW: **B5 −0.14, the best**.

### ⚠⚠ 2026-07-30: the runs were NOT 1850 GHG — a dead namelist block

**`NCMIPFIXYR: 1850` never took effect in any oifsamip run.** All 24 runscripts hand-wrote a
`&NAMECECMIP6` block. OpenIFS **48r1 renamed that namelist to `&NAMECECMIP`**
(`ece_cmip.F90:68` POSNAMs `'NAMECECMIP'` only; **`NAMECECMIP6` and `SSPNAME` appear nowhere
in the source**). The legacy block still exists in esm-tools' fort.4 *template*, so f90nml
patched it without complaint and the generated namelist looked correct on inspection — while
the model read only `NAMECECMIP` and left `NCMIPFIXYR` at its default **−1**. Since
`LFIXYEAR = NCMIPFIXYR > 0`, forcing was **transient**. Verified from the model's own output:

```
UPDRGAS: Surface greenhouse gas concentrations for YEAR/MONTH  1870/01
CO2 = 286.94 -> 287.12 -> 287.31 -> 287.49 ...   (rising through the run)
SP_SETUP: IYR = 1870, 1871, 1872, ...            (MACv2-SP aerosols advancing)
```

So every run is a **transient CMIP7 historical run at its actual calendar year**, not a
fixed-1850 PI run — Krakatoa (1883: −0.98 K Siberian JJA) and Santa María (1902) included.

**What this does NOT invalidate — the tuning results.** All 19 runs share that identical
forcing trajectory, and the perturbations live in `NAMCLDP`/`NAMVDF`/source, which *are*
read correctly. The run × year ANOVA absorbs the shared trend and the volcanoes in
`g[year]`. **B5 +0.407, B2 +0.344, B3 +0.306 stand**, as does the ±0.240 K floor.

**What it does affect — the absolute targets**, all of which are read off the control:
CO₂ drifts 287→301 ppm (≈+0.15 W/m² mid-run vs 1850, partly offset by aerosols) and the
volcanoes depress the 44-yr mean by ≈0.05 K (measured: 3 of 44 years affected). Small — but
the global TOA target is only +0.53 W/m², so ≈0.1 W/m² of spurious forcing is not negligible
*there*.

**Also unaffected: `amip_presentday`.** `historical` is *correct* for 1990–2014, so the
epoch-offset comparison is transient-vs-transient and self-consistent.

**Fix.** Not a hand-patch — the setup's own switch. `oifs.scenario: "piControl"` routes
through `configs/components/oifs/oifs48.cmip.yaml`, which writes `NCMIPFIXYR` into the live
`&NAMECECMIP` *and* sets `LCMIP_STRATAER_CMIP7=True` / `LCMIP_STRATAER_BCKGD=False` in
`NAERAD` so volcanic aerosol goes to background. `oifsamip.yaml:101` defaults
`oifs.scenario` to `"historical"`, which is why nothing was pinned. The dead block has been
removed from all 24 runscripts (no behaviour change — it was never read), and
**`amip_picontrol`** is running as a true piControl reference. `amip_pi_base` is deliberately
untouched, because all 19 lever deltas are referenced to it.

**Residual code gap.** There is **no switch to pin MACv2-SP tropospheric aerosols.**
`updtim.F90:812-818` takes the year from `UPDCAL(...NINDAT...IAN)` unconditionally;
`NCMIPFIXYR` is consulted only at `:842` for the separate `SUECSO4` sulphate path. So even
`amip_picontrol` will advance MACv2-SP 1870→1915 while its GHG, ozone and volcanic aerosol
sit at 1850. Worth reporting upstream.

**Coupled rounds 06–09 are NOT affected** — their `NCMIPFIXYR = 1850` sits in the live
`&NAMECECMIP` block, written by the proper `oifs48.cmip.yaml` machinery.

### ⚠ SUPERSEDED 2026-07-30 — the estimate below was ~2.7× too large

**Read the "Period offset, measured directly" section first.** The indirect estimate below
put the Siberian JJA period offset at ~1.12 K and the SO SW CRE offset at −1.43 W/m². The
`amip_presentday` run has since measured both: **+0.42 K** and **−0.07 W/m²**. The reasoning
below is kept because the *method* error is instructive — it applied an **observed** epoch
change to a model that does not reproduce that change — but its conclusions are withdrawn.

Measured 2026-07-29 (`scripts/analysis/era5_period_offset.sh`). The T2m reference in
`eval_round10_A.py` is `obs/era5/netcdf/T2M.nc`, which is ERA5 **1990–2014**. The AMIP runs
are **1870s observed SST with 1850 GHG**. That is ~130 years of greenhouse warming sitting
inside a number the campaign has been treating as model error and tuning against.

ERA5 Siberian JJA (55–75N, 60–180E box, no land mask), relative to the 1990–2014 reference:

| period | offset vs 1990–2014 |
|---|---:|
| 1979–1989 | **−0.670 K** |
| 1940–1969 | **−0.790 K** |

The interannual sd is 0.578 K, so the SE of the 1979–89 vs 1990–2014 difference is 0.209 K
→ `t ≈ 3.2`. The offset is real, not sampling noise. ERA5 cannot reach the 1870s, so chain
through HadCRUT5's global series (1870s → 1990–2014 = **+0.756 K**) scaled by ERA5's own
Siberian-JJA-to-global amplification over the overlapping window (**×1.48**):

```
implied Siberian JJA 1870s -> 1990-2014 offset  ~ 1.12 K
reported AMIP "bias" vs ERA5 1990-2014          = -2.16 K
residual genuine model cold bias                ~ -1.04 K
```

**The boreal target is therefore ~1.0 K, not 2.2 K.** Consequences:

1. **Tuning to close 2.2 K would manufacture a ~1 K warm bias** in the coupled PI model.
2. **The stack is smaller than feared** — at 44 yr (±0.27 K) a ~1.0 K target is ~4
   resolvable levers, not ~7.
3. **The energy target splits.** CERES EBAF here is the 07/2005–06/2015 climatology scored
   against an 1870s model. Measured in the next section: the **SO SW CRE gap is overstated**
   (~6 W/m² real, not +8.08), but the **global TOA target of ~0 is correct** and needs no
   correction, so A1b's energy result stands.

*Caveats:* the ERA5 box is unmasked while the model numbers are land-masked — land warms
faster than ocean, so the true land-only offset is if anything **larger** and this
correction is conservative. The ×1.48 amplification is extrapolated back from one window,
and HadCRUT5's 1870s global mean rests on sparse coverage.

**The clean fix is one cheap run.** The AMIP SST forcing covers 187001–201512 and CMIP6
historical GHG runs to 2014, so a present-day AMIP leg (1989–2015, transient GHG via
`NCMIPFIXYR: 0`) can be scored against ERA5 1990–2014 *and* CERES 2005–2015 with no period
mismatch at all — ~25 years, one job.

**Does NCEP2 corroborate "LPJ-GUESS was calibrated under a too-generous CRUNCEP forcing"?
On temperature, no.** NCEP2 Siberian JJA is **−0.118 K colder** than ERA5 over 1990–2014
(paired by year, sd 0.228, `t = −2.6`) — statistically real but negligible, and the *wrong
sign* for a story in which LPJG grew its trees under an over-warm forcing. This is
consistent with how CRUNCEP3 is built: temperature derives from **CRU**, bias-corrected
against station observations, while **radiation** derives from NCEP — which is exactly where
the **+21 W/m² excess against CERES** was measured. The calibration mismatch is on the
radiation axis, not the temperature axis. Caveat: NCEP2 is only a proxy for CRUNCEP3; the
direct check needs the actual CRUNCEP3 forcing, and the `.ins` files under the spin-up's
`config/lpj_guess/` are templates with placeholder paths (`c:/nc/temp.nc`), so it was not
done here. See also [[forcing-transfer-test]] — swapping CRUNCEP→AMIP forcing alone costs
−54 % Siberian TREEFPC, which is the direct evidence and does not depend on *which* variable
carries it.

### The same test on the ENERGY target: TOA survives, the SO CRE target does not

Measured 2026-07-29 (`scripts/analysis/ceres_period_offset.sh`). CERES cannot answer this —
it starts in 2000 — so ERA5's own TOA fields (178 TSR, 208 TSRC, 179 TTR) stand in, 1940–2015.

| ERA5 period | SO SW CRE | global net TOA |
|---|---:|---:|
| 1940–1969 | −59.886 | **−0.052** |
| 1970–1999 | −60.654 | +0.276 |
| 2005–2015 (the CERES window) | −61.312 | **+0.735** |

**The global TOA target needs no correction, and this is the reassuring half.** We target
~0 for a pre-industrial state, and ERA5's 1940–69 value is −0.05 — i.e. essentially
equilibrium, exactly as a near-PI climate should be. Meanwhile its 2005–2015 value of
**+0.735 W/m² reproduces the observed present-day energy imbalance (~+0.7–0.9)**, which is
a useful validation that these fields are not nonsense. So the **+0.53 W/m² model imbalance
is genuine model error**, the ~0 target is right, and **A1b's −0.51 W/m² correction stands
unqualified**.

**The SO SW CRE target does carry the flaw.** That target is taken straight from CERES
(−68.14). ERA5 says 1940–69 was **+1.425 W/m² less negative** than the CERES window
(interannual sd 1.005, SE of the difference 0.354 → `t ≈ 4.0`), and the 1870s would be
further still — call it ~+2 W/m². So the real gap for an 1870s model is **~6 W/m², not
+8.08**, and A1b's −1.89 closes **~30 % of it rather than 23 %**. The direction of the
error is unchanged; its size was overstated.

**Caveat, and it is a serious one:** ERA5's TOA radiative fluxes are *model-derived forecast*
fields, not assimilated observations. Their long-term drift partly reflects changes in the
observing system and in the forecast model, so unlike the T2m offset — which rests on an
assimilated variable — this is an order-of-magnitude bound, not a measurement. The
present-day AMIP leg is what would settle it properly.

### ⚠ Radiation vs temperature: the boreal forest is a TEMPERATURE problem

Measured 2026-07-29 (`scripts/analysis/radiation_vs_temperature_attribution.py`) from a
colleague's spin-ups, which hold CRUNCEP temperature fixed and swap only the radiation to
CERES. **AGDD5 is identical between those two arms** (775.8 Siberia, 748.0 E. Siberia,
1157.5 NH), which is the control confirming only radiation moved.

| region | radiation only (CRUNCEP→CERES) | full swap (CRUNCEP→AMIP) | radiation share |
|---|---:|---:|---:|
| NH 45N+ | −1.7 % | −25.1 % | **7 %** |
| Siberia | −8.9 % | −54.3 % | **16 %** |
| E. Siberia | −10.5 % | −66.0 % | **16 %** |

The radiation arm moves −21 W/m² while the full swap moves ~−32, so scaling linearly puts
radiation at **≤25 %**. The daily-variability variant is indistinguishable (−7.8 % vs
−8.9 %), so it is not a sub-daily distribution effect. **Temperature carries ~75–85 % of
the boreal tree collapse.**

**Consequence for the campaign: prong A's radiation work will not restore the forest.**
Closing the entire Siberian SW deficit recovers at most a quarter of the tree loss. The
forest depends on the ~1.0 K temperature bias, not on the SO/energy target. These are two
separate jobs and should stop being described as one.

**The period-mismatch also sits on the LPJG side.** CRUNCEP3 spans **1901–2015** — 20th
century — yet was used to spin up the **pre-industrial** vegetation state the coupled 1850
model restarts from. Part of the 131-GDD5 Siberian gap is therefore a legitimate PI-vs-20thC
climate difference, not model error. 131 GDD5 over a ~120-day season is ≈1.1 K, of which the
1870s→1901–2015 offset plausibly accounts for ~0.3–0.45 K, leaving ~0.6–0.8 K genuine —
consistent with the ~1.0 K derived against ERA5 1990–2014, since CRUNCEP's period is cooler
than ERA5's. See [[reference-period-offset]].

*Caveat:* the AMIP arm used the **superseded** forcing build, ~0.9 K warmer over Siberian
JJA than the corrected one, so the temperature share is if anything understated.

### ⚠ Statistical power — read before believing any boreal number below

Measured 2026-07-29 with `scripts/analysis/noise_floor.py`: a run × year ANOVA over all 19
runs, splitting each per-year value into a tuning signal, the SST-forced year excursion
shared by every run, and internal atmospheric noise. The 4-year evaluation window gives

| diagnostic | internal sd | 95 % detection threshold | runs clearing it |
|---|---:|---:|---:|
| Siberia JJA **T2m** | 0.644 K | **± 0.89 K** | **1 / 18** (only A1a, −1.15) |
| Siberia JJA **surface SW** | 4.67 W/m² | **± 6.47 W/m²** | **1 / 18** (only A1a) |
| Siberia JJA **cloud area** | 1.99 pp | **± 2.75 pp** | **0 / 18** |
| SO **SW CRE** | 0.49 W/m² | ± 0.68 W/m² | 6 / 18 |
| global **net TOA** | 0.163 W/m² | ± 0.23 W/m² | 9 / 18 |

**The energy target is well resolved; the boreal target is not resolved at all.** Every
Southern-Ocean and global-energy conclusion in this file stands. **Every boreal ranking in
the tables below is within noise** and must not be read as a result — including "B8 is the
best boreal lever", which is `t = +1.10`.

Years needed to resolve a boreal T2m signal of size Δ: **1.0 K → 3; 0.5 K → 13; 0.3 K → 36;
0.2 K → 80.** So the 4-year screen detects any lever worth ≥ 40 % of the −2.2 K target on
its own — and against the **corrected ~1.0 K target** (section above) ±0.89 K is ~90 % of
the whole target, so at 4 years we could only ever have detected a near-total fix.
**Nine boreal levers were run and none cleared it.** That is the actual boreal
result of round 10: not a ranking, but the finding that no lever tried so far is
individually large enough to measure at this run length.

**What that does *not* license.** "No single lever is big, therefore we need a bigger
lever" is one reading and it is not supported — the far more likely route to the corrected
**~1.0 K** target is a *stack* of small levers all pushed the same way (three or four at
0.3 K each would do it, not seven). These
runs cannot distinguish "each lever is ~0.3 K and they add" from "each lever is ~0". Both
are consistent with every number in this file.

**But if stacking is the plan, the measurement problem gets worse, not better.** Individual
0.3 K contributions need ~36 years each to rank, which is unaffordable across a candidate
list. Two consequences follow, and they are the practical output of this round:

1. **Lever selection must come from physical reasoning and sign confidence, not from the
   measured 4-year ranking.** Picking the top performers out of a noise-dominated table is
   the winner's curse: it selects for favourable noise, so the stack under-delivers.
   **ABB8 is exactly that experiment.** Its three components were chosen because they
   ranked well at 4 years; stacked, they gave −0.17 K.
2. **The testable unit becomes the stack, not the lever.** A stack aimed at the full ~1.0 K
   target sits right at the 4-year threshold (±0.89 K) — detectable only if it works almost
   completely — but is comfortably resolved at 44 years (±0.27 K). So build stacks from
   physics, test the stack on the extended runs, and spend attribution effort only on a
   stack that already works.

Corollary: C1/C2/E1 are **untested, not refuted**. At 4 years we cannot separate "`RLAM`
does little" from "`RLAM` does nothing", so the boundary-layer axis remains open.

### Results

**Naming.** `A1x` = Southern Ocean levers, `A2`/`Bn` = boreal-land levers, `AB` = combined.
The `Bn` labels are *atmosphere* levers within prong A — they are **not** prong B (LPJ-GUESS).
The collision is historical; read `Bn` as "boreal lever n".

| run | change | SO SW CRE | Siberia JJA T2m | global sfc flux | verdict |
|---|---|---:|---:|---:|---|
| **A1a** | `RCL_OVERLAPLIQICE` 0.65→**0.10** | **−6.51** (80 % of gap) | **−1.15 K** | −1.63 | **overshoots**; wrecks boreal |
| **A1b** | `RCL_OVERLAPLIQICE` 0.65→**0.35** | −1.89 (23 %) | −0.03 K | **−0.13** | **best so far** — energy target met (both significant); boreal unresolved, not "untouched" |
| **A2** | `RCL_KK_CLOUD_NUM_LAND` 300→150 | −0.15 | −0.18 K | +0.36 | nothing significant anywhere; boreal ~~no traction~~ unresolved |
| **expA** | `RVRSMIN(3,4)` 250→500 | −0.40 | +0.19 K | +0.41 | ~~real but far too small~~ **within noise** — sign is unconfirmed |
| **A1c** | `RDEPLIQREFDEPTH` 500→1500 m | −0.13 | +0.06 K | +0.25 | **failed** — no SO leverage; cloud-depth selectivity idea dead |
| **B1** | `DETRPEN` 0.75E-4→0.45E-4 | −0.12 | **−0.13 K** | +0.64 | **failed** — more SW yet colder; worst tropics and global RMSE |
| **B2** | `RCLDIFF_CONVI` 10→25 | −0.07 | **+0.29 K** | +0.53 | ~~best boreal so far~~ **boreal within noise**; +1.23 subpolar N Atl RMSE |
| **AB** | A1b + B2 | −1.82 | **+0.32 K** | +0.23 | **composes on energy** (SO CRE −1.82, significant); boreal within noise |
| **B3** | `RCLDIFF` 6.0E-6→1.5E-5 | **+0.54** | **+0.37 K** | **+1.52** | ~~boreal works~~ boreal within noise; **energy fails** — sfc flux +1.14, tropics +1.72 |
| **B4** | `ENTSHALP` 2.0→3.0 | −1.55 | **−0.19 K** | −1.07 | boreal within noise; **significant** SO CRE −1.55; only run to improve subpolar N Atl |
| **B5** | `RCAPDCYCL` 2.0→0.0 | −0.46 | +0.00 K | −0.58 | boreal within noise (as is every run); best global SW RMSE (13.820) |
| **B6** | `RLCRITSNOW` 2.0E-5→1.0E-5 | +0.05 | −0.03 K | +0.05 | **significant** global TOA −0.34; boreal unresolved |
| **B7** | `RVICE` 0.16→0.22 | +0.99 | **−0.67 K** | +0.06 | ~~worst boreal~~ boreal within noise; **significant** SO CRE +0.99 (wrong way) |
| **B8** | `RVLAMSK`/`RVLAMSKS(3,4)` 10→5 | +0.20 | **+0.50 K** | +0.44 | ~~best boreal of the campaign~~ **within noise** (t=+1.10); cloud −1.98 pp vs ±2.75 pp threshold |

### Deep-water-formation SW RMSE vs CERES — the priority metric
SW biases where deep water forms set the coupled ocean's initial state; errors there give
long, unpredictable coupled spin-up drift. **This ranks the runs differently from the
regional means, and it is the ranking that counts.**

**Significance-tested 2026-07-29** (`scripts/analysis/rmse_significance.py`), same run × year
ANOVA, using each year's own spatial RMSE as the replicate. Unlike the boreal means, this
metric is **well resolved** — RMSE aggregates thousands of gridcells, so internal noise
largely averages out:

| region | 95 % threshold | verdict |
|---|---:|---|
| SO 45–65S | ±0.055 | 11/18 significant — **A1b −0.16, AB −0.17, ABB8 −0.15** all real |
| subpolar N Atl | ±0.078 | 5/18 significant — see the cancellation below |
| Nordic Seas | ±0.206 | **0/18 — nothing is significant** |
| global SW | ±0.067 | 6/18 — worst B4 +0.28, best B5 −0.13 |

**The cancellation claim is CONFIRMED.** It was flagged as load-bearing and unchecked:

| subpolar N Atl | Δ | t | verdict |
|---|---:|---:|---|
| B2 alone | +0.101 | **+2.53** | **significant damage** |
| AB (A1b+B2) | +0.015 | +0.37 | within noise |
| ABB8 (A1b+B2+B8) | +0.002 | +0.05 | within noise |

B2's subpolar damage is real, and adding A1b removes it to indistinguishable from zero. The
energy and boreal targets genuinely do **not** trade against each other in that basin.

**But the Nordic Seas column below is entirely noise.** Nothing reaches ±0.206, so A1b's
−0.989 there was never established, and "A1b is best or neutral in all three deep-water
regions" holds only for two of them.

*Magnitudes differ between the two tables by construction:* the table below is the RMSE of
the 4-year **mean** field (smaller — averaging suppresses variability), while the test uses
**per-year** RMSE, the correct replicate for significance. Compare differences and signs,
not absolute levels.

| | SO 45–65S | subpolar N Atl | Nordic Seas | global SW | T2m vs ERA5 |
|---|---:|---:|---:|---:|---:|
| control | 7.193 | 4.813 | 9.231 | 14.349 | 1.585 |
| A1a | 6.179 | 5.048 | 9.617 | 14.369 | 1.578 |
| **A1b** | **5.557** | **4.819** | **8.243** | 14.138 | **1.557** |
| A2 | 6.914 | 5.052 | 9.302 | 14.275 | 1.608 |
| expA | 6.427 | 5.414 | 8.910 | **14.032** | 1.574 |
| A1c | 6.874 | 5.134 | 9.688 | 14.088 | 1.621 |
| B1 | 6.962 | 5.942 | 9.200 | 14.593 | 1.628 |
| B2 | 6.483 | **6.039** | 8.704 | 14.064 | 1.589 |
| **AB** | **5.306** | 4.885 | 9.353 | 14.041 | 1.590 |
| B3 | 6.756 | 5.888 | 8.943 | 14.375 | 1.571 |
| B4 | 6.658 | **4.440** | 8.757 | 15.303 | 1.605 |
| B5 | 6.343 | 5.387 | 8.805 | **13.820** | 1.576 |
| B6 | 6.944 | 5.209 | 9.526 | 14.067 | 1.609 |
| B7 | **7.235** | 5.425 | **9.682** | 14.471 | **1.672** |
| B8 | 7.118 | 5.464 | 9.287 | 14.221 | 1.577 |

*(Nordic Seas column: **none** of these differences is significant — see above.)*

**A1b is best or neutral in all three deep-water regions** and the only *single* lever that
leaves the subpolar North Atlantic alone. Note it beats A1a in the SO *by field RMSE*
(5.56 vs 6.18) despite A1a fixing more of the mean CRE — A1a improved the regional mean
while degrading the spatial pattern, which is the compensating-error signature showing up
in the metric built to catch it.

**The combination result looks like the important one, with one caveat.** AB gives the best
Southern Ocean of any run (5.306 vs 7.193 control) and appears to **cancel B2's
subpolar-Atlantic damage**: +1.226 for B2 alone collapses to +0.071 when A1b is added. If
real, that matters a great deal — it would mean the energy and boreal targets are not
forced to trade against each other in that basin.

**Caveat: these spatial RMSEs have not been significance-tested.** The noise floor above
was measured for regional means, not for field RMSE. Given that SO *mean* CRE is well
resolved (±0.68) the ocean RMSEs are probably sound too, but B2's subpolar signal is the
load-bearing number here and it has not been checked. Test before relying on it.

B4 is the only run that ever *improves* the subpolar North Atlantic (−0.373). Its boreal
T2m (−0.19 K) is within noise and carries no information, but its SO CRE (−1.55) and global
TOA (−1.43) are both significant, and its global SW RMSE is the worst of any run (15.303).
Filed as a possible corrective term if a future combination overshoots in that basin.

*(T2m RMSE folds in the expected PI-vs-present-day offset of ~0.5–1 K; the spread across
runs is only ±0.04 K, except B7 at +0.087. Use it to rank, not as a skill score.)*

### Round of 2026-07-28 — predicted vs actual

All eight completed; results are folded into the two tables above. The rationales are kept
because the contrast between what each lever was *supposed* to do and what it did is what
now points the campaign at the boundary layer.

| run | change | type | rationale (as written before the runs) | outcome |
|---|---|---|---|---|
| **AB** | `RCL_OVERLAPLIQICE` 0.35 **+** `RCLDIFF_CONVI` 25 | namelist | do the two halves compose? if additive: sfc flux ≈ +0.02, boreal ≈ +0.26 K | composes **on energy**; boreal +0.32 K is within noise |
| **B3** | `RCLDIFF` 6.0E-6→1.5E-5 | namelist | strongest-σ cloud knob (SPP allows ×2.83); erosion is **phase-agnostic**, so it sidesteps the mixed-phase opposition | boreal +0.37 K **within noise**; **sfc flux +1.14** — energy target lost |
| **B4** | `ENTSHALP` 2.0→3.0 | namelist | Vial/Bony: stronger shallow mixing dries the sub-cloud layer and *reduces* low cloud | boreal −0.19 K **within noise** — shallow-mixing idea untested, not dead |
| **B5** | `RCAPDCYCL` 2.0→0.0 | namelist | CAPE diurnal-cycle correction exists for **land** diurnal convection — inherent land bias | +0.00 K, but so is everything — **unresolved**, not proven inert |
| **B6** | `RLCRITSNOW` 2.0E-5→1.0E-5 | namelist | removes ice faster; targets the mid/high ice cloud that carried A1a's response | −0.03 K — **unresolved**, not proven inert |
| **B7** | `RVICE` 0.16→0.22 | namelist | same target, different route; beyond EC-Earth's 0.137–0.17 range | −0.67 K **within noise**; SO CRE +0.99 is significant and wrong-way |
| **B8** | `RVLAMSK`/`RVLAMSKS(3,4)` 10→5 | **source** | skin-layer conductivity, **vegetation-type indexed** — the only lever that cannot *directly* degrade the deep-water regions | **+0.50 K but t=+1.10 — within noise**; no mechanism can be inferred |

**B8 was called as the one to watch. That call cannot be judged from these runs.**
Its +0.50 K is `t = +1.10` and its cloud change (−1.98 pp) sits under a ±2.75 pp threshold,
so nothing about it is resolved. An earlier version of this section read the +0.50 K and the
cloud reduction as evidence that boreal cloud here is *moisture-supply-limited rather than
microphysics-limited* — inferred from a surface parameter apparently out-cutting six cloud
knobs, with A2 and B5/B6 as corroboration. **That inference is withdrawn.** It was built on
differences smaller than the noise, and A2/B5/B6 being "inert" is indistinguishable from
their being unresolved. The mechanism may still be true; these runs simply cannot speak to
it.

**What survives is one bounded statement, and only one.** Nine boreal levers spanning cloud
microphysics, convection, erosion, shallow mixing, ice fall speed and surface conductivity
each failed to move Siberian JJA T2m by the ±0.89 K that 4 years can detect. That bounds
each lever's *individual* size; it says nothing about whether they are ~0.3 K and additive
or ~0. Read it as an upper bound per lever, not as evidence that these parameters are
inert, and not as a reason to abandon incremental tuning — stacking small levers remains
the most plausible route to −2.2 K.

It also leaves ECMWF's own diagnosis — excessive turbulent mixing in cloudy boundary layers,
explicitly *not* cloud microphysics — as the leading untested hypothesis rather than a
tested one. C1/C2 were the first shot at it and returned no resolvable signal, which at
4 years means untested, not refuted.

### Round of 2026-07-29 — all four completed, none resolved on boreal

Reference for **C1/C2/E1 is `amip_B8_lamsk5`**, not the control: the current `install/lib`
has B8 compiled in, so every run launched on it inherits `RVLAMSK/S(3,4)=5`. Their
namelists are otherwise byte-identical to the control's (`RVICE: 0.16` only), so each is a
clean single-lever increment on a known baseline.

| run | change | type | binary | rationale |
|---|---|---|---|---|
| **ABB8** | A1b 0.35 + B2 25 + B8 | namelist on B8 | `e41a8d4acdb3` | do the three compose? → **energy yes** (SO CRE −1.95, TOA −0.27, both significant); boreal −0.17 K, within noise |
| **C1** | `RLAM` 150→**75** m | namelist (`NAMVDF`) | `e41a8d4acdb3` | BL mixing, new axis → −0.21 K, **within noise: untested, not refuted** |
| **C2** | `RLAM` 150→**40** m | namelist (`NAMVDF`) | `e41a8d4acdb3` | aggressive → +0.39 K, also within noise; C1/C2 non-monotonicity is noise, not physics |
| **E1** | `RVLAMSK`/`RVLAMSKS(3,4)` 5→**2.5** | **source** | `4b1f678bc051` | linear or saturating? → **unanswerable**: +0.06 K, and B8's +0.50 K was itself unresolved |

**The C group is the important new idea.** `RLAM` is the asymptotic mixing length used
**only in the statically-unstable branch** of `vdfexcu` (`vdfexcu.F90:397`, `ZKLENT=PLAM`)
— the daytime convective land boundary layer and nothing else. `suvdf.F90:77` sets the
150 m default and `:108-109` then reads `NAMVDF` over it, so it is namelist-only, needs no
rebuild, and can run in parallel with source-edit experiments.

It matters because it is the one lever aimed at ECMWF's *own* published diagnosis of the
IFS summer land cold bias: **excessive turbulent mixing in cloudy boundary layers**, with
the explicit statement that cloud microphysics is *not* the main driver. Everything in
rounds A and B was a cloud or convection knob; the whole B series returned at most +0.5 K,
which is consistent with ECMWF being right and us having spent the campaign on the wrong
axis. Less mixing → shallower BL → surface sensible heating retained in a thinner layer →
warmer T2m, *and* less moisture lofted → less boreal cloud. Both signs favour us.

Caveat: `RLAM` is global and applies over ocean too, so the deep-water SW RMSE must be
checked, not assumed.

**E1 was posed as a decision run — linear or saturating? — and the honest answer is that
the question was ill-posed.** It rested on B8's +0.502 K being a measured quantity, and it
was not (`t = +1.10`). E1 returned +0.06 K, equally unresolved. Nothing about the skin-
conductivity axis can be concluded from B8 and E1 together, in either direction.

The design error worth remembering: **B8 was promoted to "best of campaign" and then built
on — twice, in ABB8 and E1 — before its uncertainty was ever estimated.** Three of the four
runs in this round were spent following up a result that was inside the noise. The noise
floor cost one post-processing script and no node-hours, and should have been measured
after the first round, not the third.

**Serialisation cost.** Source-edit experiments share one model tree, so only one can be
in flight at a time (edit → `recomp` → verify md5 → submit → wait for staging → repeat).
Namelist experiments have no such constraint. If the E/vegetation axis proves worth
several more runs, copying the model tree would let them run in parallel — worth raising
before doing it, given the stale-object incident below.

### BUILD PROCEDURE — non-negotiable for source edits
`esm_master comp-` **silently reused a stale object** for `susveg_mod.F90`: the object md5
was byte-identical with `RVLAMSK=5` and `RVLAMSK=10`. Almost certainly because
`source/ifs_sp` is a symlink to `ifs-source`, so CMake's dependency check follows the link,
not the edited file. This cost five wrongly-killed runs and a wrong contamination call.

1. edit source
2. **`esm_master recomp-oifsamip-cy48/oifs-48r1`** from `model_codes` (`recomp` = conf +
   clean + comp; the `clean` is what forces it). `comp-` is fine for namelist-only work.
3. **verify the object or library md5 actually changed** — not the source, the *binary*
4. only then submit
5. **do not rebuild while a run is queued or starting**: `esm_tools` copies `install/lib/*`
   into `work/lib/oifs` at **job start**, not at submit time. Wait until the run's
   `work/lib/oifs/libsurf.SP.so` exists, then rebuild for the next experiment.

### What has been learned

**1. `RCL_OVERLAPLIQICE` = 0.35 essentially solves the energy target.** Global surface
flux +0.383 → **−0.129 W/m²** and net TOA +0.533 → +0.027, at a cost of only −0.26 in the
tropics and no change in NH−SH albedo. One namelist parameter, no rebuild.

**2. EC-Earth4's 0.1 is not transferable to our configuration.** It fixes 80 % of the SO
CRE error but drives global TOA to −1.48 and cools Siberian JJA by 1.15 K. EC-Earth4
presumably runs it alongside compensating tuning we do not have.

**3. The A1/A2 separability argument was WRONG, and conditional on amplitude.** The
prediction was that the WBF terms could not touch boreal summer cloud because
boreal-land BL cloud is warm (>268 K) and never enters the mixed-phase window. Measured:
it **fails badly at 0.1 (−1.15 K, significant)**. At 0.35 the response is −0.03 K, which is
unresolved — consistent with the prediction but not a confirmation of it. **The error was ignoring the
mid-level mixed-phase cloud above the warm boundary layer.** Vertical decomposition of
A1a over Siberian land in JJA:

| | control | A1a change |
|---|---:|---:|
| low cloud | 52.3 % | +1.06 |
| **mid cloud** | 39.8 % | **+3.50** |
| high cloud | 46.4 % | +0.52 |
| **column liquid water** | 72.3 g/m² | **+14.6 (+20 %)** |
| column ice water | 23.9 g/m² | −4.2 |

Lowering the overlap converted ice→liquid in mid-level cloud; liquid is far more
reflective per unit mass, hence −9.8 W/m² surface SW and the cooling.

*Caveat on this table:* the ±2.75 pp cloud-area detection threshold applies here too, so the
individual layer changes are at best marginal on their own. What the mechanism has to
explain — A1a's −1.15 K and −9.75 W/m² — is significant, and the +14.6 g/m² column-liquid
change is large; but treat the layer-by-layer split as indicative rather than measured.

**4. This sets up a direct opposition through every mixed-phase knob.** The Southern
Ocean needs *more* supercooled liquid (brighter); boreal mid-level cloud needs *less*
(dimmer). Both are the same process. **No mixed-phase parameter can fix both**, which is
why A1 helped the SO and hurt the boreal. Any boreal fix must act through a
non-phase-partitioning channel.

**5. The boreal cold bias is the hard problem, and it is also the badly-measured one.**
expA +0.19 K and A2 −0.18 K are both **inside the ±0.89 K noise floor**, so neither the
sizes nor the *signs* are established — an earlier version of this entry called expA
"real" and A2 "wrong sign", and neither is supported. What does stand: A1 makes the boreal
worse at any strength that meaningfully fixes the SO (A1a −1.15 K **is** significant), and
the energy target fell to a single parameter.

**6. Warm-rain removal did not detectably control boreal cloud.** A2 sped land
autoconversion 3.5× and moved JJA cloud area by −0.08 pp — but the detection threshold on
that diagnostic is **±2.75 pp**, so this is an upper bound, not a measurement. The
conclusion once drawn from it — that boreal summer cloud is moisture-supply-limited rather
than rainout-limited — is **not supported by this run**. It remains a plausible hypothesis
with no evidence behind it either way.

### Standing caution
A1b reaches near-zero global flux while leaving **77 % of the SO CRE error and the entire
boreal error in place**. That is precisely the compensating-error configuration
Schuddeboom & McDonald (2021) warn about — CMIP6 models with the *smallest* SO radiation
bias carry the *largest* compensating errors. The global number looking right is not
evidence the cloud field is right. Tuning to 0.30–0.35 to hit −0.16 exactly would be
tuning to the global number, which is the failure mode the literature names.

---

Source tree: `/work/ab0246/a270092/model_codes/oifsamip-cy48/oifs-48r1/ifs-source/`
(branch `movcav-landice+co2-concdriven` @ `f3ccacb`, = esm_tools version `48r1v5`).

**What we currently override:** `RVICE = 0.16`, `GGAUSSB = -0.5`, `ENTSTPC3 = 1`.
Everything else runs at OpenIFS defaults — i.e. **our configuration is not the
EC-Earth-tuned one**. Rounds 06–09 tuned ocean and sea-ice parameters on top of a
largely untuned atmosphere.

---

## 1. The headline finding: `RCL_OVERLAPLIQICE`

| | value |
|---|---|
| OpenIFS 48r1 default (`arpifs/phys_ec/sucldp.F90:296`) | **0.65** |
| EC-Earth4 `tuning_v0.3.yml` — their *only* OIFS tuning entry | **0.1** |
| ours | 0.65 (unset) |

Assumed sub-grid overlap fraction of supercooled liquid with ice in the
**Wegener–Bergeron–Findeisen deposition** term
(`cloudsc.F90:2499 → 2540`, and `:2603` for `IDEPICE=2`):

```
ZDEPOS = MAX(ZOVERLAP_LIQICE * ZA * (ZINEW - ZICE0), 0)
```

Lowering it exposes less liquid to growing ice crystals → supercooled liquid survives →
brighter cloud. This is the mechanism the Southern Ocean literature blames for the SW
bias, and EC-Earth4 converged on it independently.

ECMWF also wrote a **regime-selective version and left it commented out**
(`cloudsc.F90:2494-2504`, mirrored at `:2605-2608`) — verified directly:

```fortran
! Reduce in shallow convection because assume SLW in active
! updraught is less overlapped with ice in less active part
ZOVERLAP_LIQICE = RCL_OVERLAPLIQICE
!IF (KTYPE(JL) > 0 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
!  ZOVERLAP_LIQICE = 0.1_JPRB
!ENDIF
```

Same 0.1 EC-Earth4 chose globally. Re-enabling it, ideally regated on `PLSM < 0.5` or on
estimated inversion strength `PEIS` rather than `KTYPE`, is the cheapest route to a
genuinely Southern-Ocean-only version.

---

## 2. Mixed-phase partitioning in OIFS — not what the literature assumes

Important structural correction. The usual "models put too much ice in at condensation"
story **does not apply to OIFS 48r1**:

- **Condensate phase is a hard step at `RTHOMO` = RTT − 38 = 235.15 K**
  (`cloud_satadj.F90:628-632`, `:669-673`, `:762-770`). Stratiform condensate is
  **100 % liquid down to −38 °C**.
- The `FOEALFA` quadratic ramp between `RTICE` (=RTT−23=250.16 K) and `RTWAT` (=RTT=273.16 K)
  — `fcttre.func.h:83-84`, constants in `suphec.F90:193,196` — only builds the mixed-phase
  *saturation vapour pressure* `ZQSMIX`. It does **not** partition condensate.
- Splitting of *existing* condensate uses the prognostic ratio `qL/(qL+qI)`
  (`cloudsc.F90:1266-1267`), not temperature.
- Convectively detrained condensate uses a *third* ramp `FOEALFCU` with
  `RTICECU`, which `LMFGLAC=.TRUE.` (default) overwrites to RTT−38 (`sucumf.F90:246`).
  Shallow convection detrains **100 % liquid** regardless of T (`cuascn.F90:820`).

**Consequence:** excess ice comes entirely from the **conversion** terms — WBF deposition,
riming, snow autoconversion — not from the condensation partition. `RTICE`/`RTWAT` are
hardcoded and, in any case, the wrong target.

---

## 3. NAMCLDP — namelist-settable cloud parameters

Defaults in `arpifs/phys_ec/sucldp.F90`, all set before the namelist read (line 627), so
all are overridable. Sign = effect of **increasing** the parameter on cloud reflectivity.

| Parameter | Default | file:line | Effect | Sign |
|---|---|---|---|---|
| `RCL_OVERLAPLIQICE` | **0.65** | :296 | WBF liquid/ice overlap — see §1 | **−** |
| `RCLDIFF` | 6.0E-6 s⁻¹ | :247 | Turbulent erosion at cloud edge | − |
| `RCLDIFF_CONVI` | 10.0 | :248 | Multiplier on `RCLDIFF` at convective points; **the only regime-gated NAMCLDP knob** (see §6) | − |
| `RVICE` | 0.13 m/s *(we use 0.16)* | :273 | Cloud-ice terminal fall speed | − |
| `RLCRITSNOW` | 2.0E-5 kg/kg | :267 | Critical in-cloud ice for ice→snow autoconversion | + |
| `RSNOWLIN2` | 0.030 K⁻¹ | :266 | T-dependence of ice→snow rate (`0.025` = Lin et al. 83) | + (cold only) |
| `RDEPLIQREFRATE` | 0.5 | :294 | Fraction of deposition rate in cloud-top layer | − |
| `RDEPLIQREFDEPTH` | 500.0 m | :295 | Depth of the supercooled-liquid cloud-top layer | + |
| `RCL_EFFRIME` | 1.0 | :301 | Riming efficiency — comment says physical range is **< 1**, so it sits at its maximum | − |
| `RCL_INHOMOGAUT` | 1.5 | :257 | Sub-grid inhomogeneity enhancement on KK-2000 **autoconversion** | − |
| `RCL_INHOMOGACC` | 3.0 | :258 | Same for **accretion** | − |
| `RKOOPTAU` | 10800 s | :330 | Ice-supersaturation removal timescale | − |
| `RTAUMEL` | 7200 s | :325 | Snow/ice melting relaxation | ~neutral |

**`RCL_INHOMOGAUT` = 1.5 and `RCL_INHOMOGACC` = 3.0 are the "warm-rain enhancement
factors" the literature flags** (1.5× autoconversion, 3× accretion on Khairoutdinov–Kogan
2000) — reported worth ~35 W/m² downward SW and 40–60 g/m² cloud water in an IFS-family
model. Both are namelist-settable here.

### Dead knobs — do not spend time on these
| Parameter | Why inert |
|---|---|
| `RKCONV` | Only used under `IWARMRAIN==1`; `IWARMRAIN=3` is hardwired (`cloudsc.F90:752`) |
| `RCLCRIT_SEA`/`RCLCRIT_LAND` | Computed but unused in the `IWARMRAIN==3` branch |
| `RCCNOM`/`RCCNSS`/`RCCNSU` | Gated off by `NAERCLD=0` |
| `RTAU_CLD_TLAD` | TL/AD only |
| `RSWINHF`/`RLWINHF` | Inert unless `NINHOM=1` |
| `RCLOUD_SEPARATION_SCALE_*` | SPARTACUS only; inert at `NSWSOLVER=0` |

*(`RCLD` and `RTINT` do not exist in this tree.)*

---

## 4. NAMCUMF — convection (`sucumf.F90`)

| Parameter | Default | file:line | Effect |
|---|---|---|---|
| `ENTRORG` | 1.75E-3 m⁻¹ | :122 | Organized entrainment, positively buoyant convection |
| `ENTSHALP` | 2.0 | :125 | Shallow entrainment = `ENTSHALP × ENTRORG` |
| `DETRPEN` | 0.75E-4 m⁻¹ | :115 | Detrainment for penetrative convection → sets `PLUDE` fed to `cloudsc` |
| `RPRCON` | 1.4E-3 | :156 | Updraught water → precipitation conversion |
| `ENTRDD` | 3.0E-4 m⁻¹ | :136 | Downdraught entrainment |
| `RMFDEPS` | 0.30 | :146 | Fractional downdraught mass flux at LFS |
| `ENTSTPC3` | 3.0 *(we use 1.0)* | :131 | Extra entrainment — **only** for the PBL inversion height |
| `RDEPTHS` | 2.E4 Pa | :151 | Max shallow cloud depth |
| `RTAUA` | 1.0 | :174 | CAPE-closure timescale multiplier |
| `RHEBC` | 0.92 | :171 | Critical RH below cloud for downdraught evaporation |
| `RCAPDCYCL` | 2.0 | :208 | CAPE diurnal-cycle correction |
| `LMFGLAC` | .TRUE. | — | `.FALSE.` restores `RTICECU`=RTT−23, steepening the detrained-phase ramp |

Hardcoded (not in NAMCUMF): `RMFCMIN`, `RCUCOV`, `RCVRFACTOR`, `LMFMID`, `LMFUVDIS`,
`LMFWSTAR`, `LMFSMOOTH`, `LMFPROFP`.

---

## 5. NAERAD — radiation (`suecrad.F90`). No `NAMECETUNING`/`NAMRAD` in this tree.

| Parameter | Default | file:line | Effect |
|---|---|---|---|
| `RCLOUD_FRAC_STD` | 1.0 | :776 | In-cloud water FSD for McICA |
| `RMINICE` | 60.0 µm | :744 | Min ice effective diameter; with `NMINICE=1`, `Dmin = 20 + (RMINICE−20)·cos(lat)` |
| `NMINICE` | 1 | :743 | 0 = constant `Dmin`; 1 = latitude-varying |
| `NRADIP` / `NRADLP` | 3 / 2 | :737/:738 | Ice (Sun 2001) / liquid (Martin 1994) effective radius |
| `RCCNSEA` / `RCCNLND` | 50 / 900 cm⁻³ | :861/:860 | Marine/land CCN — **inert** unless `LCCNO`/`LCCNL` set `.FALSE.` |
| `NCLOUDOVERLAP` | 3 | :690 | Exponential-random overlap |

---

## 6. Land/sea and regime selectivity — the honest answer

**No NAMCLDP parameter is land/sea selective.** Every one is a global scalar. Ranked by
what selectivity is achievable:

1. **The mixed-phase temperature window is the best de-facto selectivity.** Deposition,
   riming and snow-autoconversion terms only act where supercooled liquid coexists with
   ice (`RTHOMO ≤ T < RTT−5 K`, `qL > RLMIN`, `cloudsc.F90:2512`). **Boreal-land summer
   BL cloud is warm (T > 268 K) and sees essentially none of them.** So an SO fix via
   `RCL_OVERLAPLIQICE` / `RDEPLIQREF*` / `RCL_EFFRIME` should leave boreal JJA alone —
   a falsifiable prediction worth checking in every A1 run. Collateral falls on
   mid-latitude and Arctic *winter* mixed-phase cloud, boreal winter/spring, and the
   mid-troposphere.
2. **`RCL_KK_CLOUD_NUM_SEA` = 50 cm⁻³ (`sucldp.F90:417`)** vs `..._LAND` = 300 (`:418`) —
   genuinely ocean-only, selected by `IF (PLSM > 0.5)` at `cloudsc.F90:2140-2144`, entering
   KK autoconversion as `N^(−1.79)`. **Hardcoded** — needs a source edit. Raising 50→80
   slows marine warm-rain autoconversion ~2.3× with zero land effect. No radiative effect
   through droplet size (radiation gets N_d independently from Martin + climatological CCN).
3. **`RCLDIFF_CONVI`** — the only regime-gated NAMCLDP knob. Under `ITURBEROSION=3`:
   `20×` if `KTYPE≥2` **and** `PEIS < 10 K` (weak inversion), `2×` if `KTYPE==1`
   (`cloudsc.F90:1779-1783`). Hits weak-inversion cumulus-topped regimes — including SO
   cold-air outbreaks — harder than capped stratocumulus. Still touches continental
   shallow convection.
4. **`RCCNSEA` with `LCCNO=.FALSE.`** — namelist-settable and ocean-only, but replaces a
   spatially varying climatological CCN field with one global number. Structural, not a nudge.

---

## 7. Hardwired scheme selectors (`cloudsc.F90`) — recompile to change

`IWARMRAIN=3` (:752), `IRAINACC=1` (:759), `ISUBLSNOW=1` (:771), `IDEPSNOW=1` (:777),
`IDEPICE=1` Rotstayn 2001 (:783), `ISUBLICE=0` (:789), `IVARFALL=1` (:810),
`ITURBEROSION=3` (:817), `IFTLIQICE=1` (:3850).

---

## 8. ECMWF's own uncertainty magnitudes — principled tuning ranges

From the SPP (stochastically perturbed parameterizations) log-normal 1σ in
`arpifs/module/spp_def_mod.F90:246-266`. `ln1=.true.` → multiplicative, 1σ range = ×/÷ exp(σ).

| Parameter | σ | 1σ range |
|---|---|---|
| `RCLDIFF` | 1.04 | ×/÷ 2.83 |
| `RLCRITSNOW` | 0.78 | ×/÷ 2.18 |
| `RTAUA` | 0.78 | ×/÷ 2.18 |
| liquid/ice effective radius | 0.78 | ×/÷ 2.18 |
| `RAINEVAP`, `SNOWSUBLIM` | 0.65 | ×/÷ 1.92 |
| `RPRCON` | 0.52 | ×/÷ 1.68 |
| `ENTRORG`, `ENTSHALP`, `ENTSTPC1`, `DETRPEN` | 0.39 | ×/÷ 1.48 |
| `RCL_INHOMOGAUT` / `ACC` | 0.30 | ×/÷ 1.35 |
| `RAMID` | 0.13 | ×/÷ 1.14 |

Useful as defensible bounds when proposing a change.

---

## 9. EC-Earth's operational tuning sets (SMHI GitLab)

Auth: `~/.smhi_api_token` (works, user `jan.streffing`); `~/.smhi_token` is dead (401).
Group `ec-earth` = id 378. **Group-scope blob search is disabled** (no Elasticsearch);
**project-scope** search works: `/api/v4/projects/:id/search?scope=blobs&search=TERM`.

- `ec-earth/ecearth3` = id **2370** → `runtime/classic/ctrl/ifs-tuning-parameters-*.sh`
- `struthers/ecearth4` = id **2363** → `scripts/runtime/templates/tuning*.yml`
- `ec-earth/ecearth4-tuning` = id **1657** — *dedicated tuning project, not yet explored*
- `ec-earth/ec-earth-documents` = id 2411, `ec-earth/ec-earth-cmip7-experiments` = id 2210

### EC-Earth3 IFS tuning by resolution

| | T159L62 | T255L91 | T255L91-AerChem | T511L91 *(untuned)* | ours |
|---|---|---|---|---|---|
| `RPRCON` | 1.24E-3 | 1.34E-3 | 1.34E-3 | 1.2E-3 | default 1.4E-3 |
| `RVICE` | 0.17 | 0.137 | 0.137 | 0.15 | **0.16** |
| `RLCRITSNOW` | 3.50E-5 | 4.0E-5 | 4.0E-5 | 4.0E-5 | default 2.0E-5 |
| `RSNOWLIN2` | 0.035 | 0.035 | 0.03 | 0.035 | default 0.030 |
| `ENTRORG` | 1.80E-4 | 1.70E-4 | 1.75E-4 | 1.80E-4 | default 1.75E-3 |
| `DETRPEN` | 0.85E-4 | 0.75E-4 | 0.75E-4 | 0.85E-4 | default 0.75E-4 |
| `ENTRDD` | 3.0E-4 | 3.0E-4 | 3.0E-4 | 3.0E-4 | default |
| `RMFDEPS` | 0.3 | 0.3 | 0.3 | 0.3 | default |
| `RCLDIFF` | 3.E-6 | 3.E-6 | 3.E-6 | 3.E-6 | default 6.0E-6 |
| `RCLDIFFC` | 5.0 | 5.0 | 5.0 | 5.0 | default 10.0 |
| `RLCRIT_UPHYS` | 0.92e-5 | 0.875e-5 | 0.875e-5 | 0.935e-5 | — |

Notes: tuning is redone **per resolution**, not transferred. **TCO95 (~100 km) sits between
T159 (~125 km) and T255 (~80 km)**, so T159 is the closer analogue. The T511 file carries an
explicit `!!! THE T511L91 HAS NOT BEEN TUNED !!!` banner. Headers record that the tuning is
tied to a **specific LPJ-GUESS vegetation dataset** (v16 for T255, v29/ERA20C for T159) and
to `NCLOUDACT=2, NAERCLD=9` — so their numbers are **not liftable wholesale**. Provenance
tickets: `dev.ec-earth.org/issues/548` (T159), `#449` (T255).

Note EC-Earth3 halves `RCLDIFF` (3E-6 vs 6E-6 default) and halves `RCLDIFFC` (5 vs 10).

### EC-Earth4 (`tuning_v0.3.yml`) — OpenIFS 48r1, our cycle
```yaml
oifs:  tuning: namcldp: RCL_OVERLAPLIQICE: 0.1
nemo:  tuning: oce: namzdf_tke: {nn_etau: 0, rn_lc: 0.2}
```
`tuning-example.yml` additionally documents RPRCON, ENTRORG, DETRPEN, ENTRDD, RMFDEPS,
RVICE (0.13), RLCRITSNOW (0.3E-4), RSNOWLIN2 (0.3E-01), RCLDIFF (0.3E-05), RCLDIFF_CONVI (7.0).

### RPRCON: the prior AWI history (recovered 2026-08-10, project_management #87 / #95)

Written down here so `preflight.py RPRCON` surfaces it — the parameter reads as "never
set in any run", which is true of **this** campaign and false of AWI-CM3 as a whole.

| where | cycle | value | outcome |
|---|---|---|---|
| TCO319 DART, #87 (Semmler 2022) | 43r3 | 1.4E-3 → **0.7E-3** | adopted, carried through DARTC/D/E/F/J/K |
| TCO319 DART, #87 from DARTL on | 43r3 | **1.0E-3** | relaxed from 0.7E-3 |
| TCO95L91-CORE2, #95 (Streffing 2022) | 43r3 | **1.0E-3** | launched, **result never posted**, not in the AWI-CM3.1 final config |
| TCO95L91 AMIP, this campaign | **48r1** | 1.4E-3 (default) | **never set** |

What #87 measured at TCO319, quoting Semmler directly:

* *"Positive impact of BV smoothing and reduction of RPRCON visible in the tropics"* (precip vs GPCP)
* *"Slightly more clouds through reduction of RPRCON, **also over Southern Ocean**. However, no
  positive impact on the subtropical stratocumulus decks."*
* *"RPRCON reduction helps in large areas to reduce the low **liquid water path** bias"*
* *"...larger areas with positive net longwave radiation biases ... consistent with the changes
  in liquid water path"*
* Energy, DARTZ → DARTC: TOA **−0.29 → −0.71**. *"Less efficient conversion of liquid water into
  rain leads to less energy absorption of the ocean, probably because more cloud water in the
  atmosphere."*

**Two things this kills and one it opens.**

1. It is **not** a deep-convection-selective lever. RPRCON acts on the whole mass-flux scheme,
   shallow and mid-level included, and #87 measured it changing Southern Ocean cloud. The
   "keyed on deep convection so it cannot reach the SO" argument applies to `DETRPEN`, not to
   this. Do not propose RPRCON as the tropical-LW lever on that reasoning.
2. It works through **liquid** water and **removes** energy globally (SW reflection beats the
   LW gain), so it is not a tropical *warming* lever either.
3. But it demonstrably **adds Southern Ocean cloud**, and the SO cloud-*amount* deficit
   (~6 pp, two thirds of the band error) is the one term in this campaign with no identified
   lever at all. That is the hypothesis worth testing here.

**Cycle caveat, and it is decisive.** All of the above is 43r3. Project issue #170 is precisely
the finding that 48r1 has substantially different cloud cover from 43r3, so none of these
numbers transfer — they establish direction and a defensible value range (0.7–1.0E-3), not an
expected magnitude.

---

## 10. Tuning runs

All AMIP, TCO95L91, 1850 GHG, observed SST 1870s, `/work/bb1469/a270092/runtime/oifsamip-cy48/`.
Runscripts in `~/esm_tools/runscripts/oifsamip/`. All evaluated on **1872–75** against
`amip_pi_base` (discard 1870–71 for deep-soil spin-up). 6 yr ≈ 1.3 h on 14 nodes.

| run | lever | change | how | status | result |
|---|---|---|---|---|---|
| `amip_pi_base` | — reference | OIFS defaults + `RVICE=0.16` | — | **done** | net TOA **+0.67**, sfc +0.52; Siberia JJA **−2.16 K** vs CRUNCEP3, **−11 W/m²** SW vs CERES; SO cloud **−6.6 pp**, SW CRE **+7.8** |
| `amip_expA_rvrsmin500` | `RVRSMIN(3,4)` | 250 → 500 | source | **done** | **+0.19 K** Siberia JJA; SH −3.9 / LH +2.7 W/m² (Bowen 0.43→0.52); `lcc` −0.018; ΔTOA **+0.012**. Right sign, <10 % of target — mechanism confirmed, not a lever |
| `amip_A1_overlap01` | `RCL_OVERLAPLIQICE` | 0.65 → **0.1** | namelist | running | — |
| `amip_A1_overlap035` | `RCL_OVERLAPLIQICE` | 0.65 → 0.35 | namelist | running | — |
| `amip_A2_kknumland150` | `RCL_KK_CLOUD_NUM_LAND` | 300 → **150** cm⁻³ | source | pending | — |

**Falsifiable predictions** (the real test of the two-pronged design):
A1 should leave **boreal JJA unchanged** (mixed-phase window only; boreal summer cloud is
warm). A2 should leave the **Southern Ocean unchanged** (`PLSM>0.5` branch only). If either
bleeds, the separability argument fails.

**Gates for every A-run:** SO 45–65 °S SW CRE + cloud area vs CERES (**radiation only** —
prescribed SST makes an SO temperature response impossible); global net TOA *and* surface
flux (state the convention); **tropics** (50 % of the globe, largest single contribution,
and what the 06O/06T/06V levers wrecked); **NH−SH albedo** (A1 and A2 push it the *same*
way — they do not cancel on this metric); CMPI by pattern.

### Build discipline (learned the hard way)
- Use `esm_master comp-oifsamip-cy48/oifs-48r1` **from `/work/ab0246/a270092/model_codes`**;
  escalate to `recomp-…` (= conf + clean + comp, full rebuild) only if that is not enough.
  **Never** call `comp-oifs-48r1_script.sh` directly — doing so left `build/` and `install/`
  inconsistent, and runs stage `install/lib/*`, so the wrong library would have been used.
- **Judge a build only after the process has exited.** Check the log for the final
  `cp  oifs-48r1/install/bin/OpenIFS  bin` line. A mid-link timestamp comparison falsely
  reads "stale" and cost one needless full rebuild.
- Verify a source-constant change in the **object file**, not the shared library or the
  driver binary: e.g. `sucldp.F90.o` → `150.0` ×2, `300.0` ×0. Whole-library counts are
  noise, and `bin/OpenIFS` is a thin driver that does not contain the physics at all.
- Runs stage their own copy of `install/lib/*` into `run_*/work/lib/oifs/` at submit time,
  so a rebuild does **not** disturb already-submitted jobs (verified: 36 `.so` staged).

### Coupled tuning rounds 06–09
Not repeated here — see report Part I (App. A parameter catalogue) and §Round 09.
Headline: those rounds tuned **ocean and sea-ice** parameters on top of an atmosphere left
at OpenIFS defaults, which is why §9 above matters.

---

## 11. Literature leads (see report Part IV for the full argument)

- **Southern Ocean too-clear is canonical**; our −6.6 pp / +7.8 W/m² sits mid-range of
  reported values. **NH-land too-cloudy runs *against* the CMIP consensus** (models are
  usually too clear/too warm there) — so it is family-specific.
- **ECMWF's own diagnosis conflicts with ours**: they attribute the IFS summer land cold
  bias (1–2 K, Europe) to *excessive* turbulent mixing in cloudy BLs and say "cloudiness
  errors do not appear to play a major role."
- **Hemispheric albedo symmetry**: brightening the SO and dimming NH land push NH−SH the
  **same** way — they add. Track NH−SH albedo as an explicit constraint.
- **Compensating errors**: Schuddeboom & McDonald 2021 find CMIP6 models with the
  *smallest* SO radiation bias carry the *largest* compensating errors. Our global CRE
  matching CERES to 0.08 W/m² is a warning sign, not a health check.
- **Kay et al. 2016**: fixing SO SW raised ECS 4.1 → 5.6 K and needed a deliberate tropical
  dimming to hold global balance.
- **Varma et al. 2020**: recovered ~4 W/m² annual (8 seasonal) SO SW CRE via ice capacitance
  1.0→0.5 and ice-nucleation T −10→−20/−40 °C — but TOA LW loss exceeded the SW gain and
  the tropics degraded.
- **Bodas-Salcedo**: in HadGEM3 the SO SW bias is mostly cloud *cover*, not albedo —
  consistent with our −6.6 pp area signature, and an argument against phase partitioning
  being the highest-yield fix.
