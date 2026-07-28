# OpenIFS 48r1 tuning knobs — reference

Compiled 2026-07-28 from three parallel sweeps: the OIFS source tree, the SMHI GitLab
EC-Earth repositories, and the published literature. Kept as a reference because we will
likely come back here.

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
