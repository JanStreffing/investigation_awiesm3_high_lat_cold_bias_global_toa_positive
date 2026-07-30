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
