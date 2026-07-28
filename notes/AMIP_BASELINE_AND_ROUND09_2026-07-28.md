# AMIP atmosphere-only baseline, round 09, and the snow-enthalpy convention (2026-07-28)

Recovered and consolidated from the working session of 2026-07-27/28. The analysis
scripts behind every number here were run from an ephemeral scratchpad that has since
been wiped; they have been recovered from the session transcript into
`scripts/analysis/` (`budget.py`, `jja.py`, `seasonal.py`, `seasonal2.py`, `plot2.py`,
`plot_campaign.py`, `plot_toa.py`, `iav.py`).

---

## 1. New runs

| run | config | path | state |
|---|---|---|---|
| `amip_pi_base` | OpenIFS-AMIP TCO95, 1850 GHG, observed SST 1870–79, **no ocean, no LPJG** | `/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_base` | 10 yr, complete |
| `Tuning_test_09C_06V_CRUNCEPinit_newSeaIce` | 06V (06T + `ENTSTPC3=1`) on FESOM `awiesm3-implicit-ice-surftemp` + OIFS `movcav-landice+co2-concdriven` | `/work/bb1469/a270092/runtime/awiesm3-v3.4/` | 30 yr, complete |

Runscripts: `~/esm_tools/runscripts/oifsamip/oifsamip-cy48-levante-TCO95L91_amipbase.yaml`
and `~/esm_tools/runscripts/awiesm3/develop/cmip7/100Y_..._3rd_round_09C_06V_CRUNCEPinit_newSeaIce.yaml`.

Two configuration points that matter for comparability:

* The oifsamip OIFS was moved from `local_combined_fixes` @ `1303a90` to
  `movcav-landice+co2-concdriven` @ `f3ccacb` — the same branch the coupled campaign
  uses (esm_tools version `48r1v5`). The branches differ in real atmospheric physics
  (`suecozv`, `suecrad`, `suphec`, `srfi_mod`, `surfece`, `voskin_mod`), so the old
  AMIP build was **not** a valid attribution reference for the coupled model.
* `TMPRALBSEAD = 0.075` was **removed** from the AMIP runscript. The coupled runs never
  set it and therefore use the source default `RALBSEAD = 0.06`
  (`susrad_mod.F90:107`); OASIS exchanges only the *ice* albedo (`sia_feom` →
  `A_Ice_albedo`), so OpenIFS computes open-water albedo itself in coupled mode too.
  The 0.015 offset is worth order 1 W m⁻² of reflected SW — comparable to the
  imbalance being attributed.

---

## 2. Snow enthalpy: what `part2` includes, and what these numbers mean

`part2_rad_balance.py:255-268`:

```python
surface = ssr + str + sshf + slhf − sf·ρ·Lf     # snow enthalpy enters HERE, only here
toa     = tsr + ttr                              # no snow term
rad_balance = toa − surface
```

`sf` is in metres water-equivalent in these runs, so the `ρ·Lf = 3.3355e8` branch
applies after dividing by the 3600 s accumulation period.

**Every net-TOA number quoted in this campaign is `(tsr+ttr)/accumulation_period`,
i.e. exactly part2's TOA total — no snow term, in mine or part2's.** That is correct
rather than an omission: snowfall enthalpy is an internal atmosphere↔surface transfer
and cancels at the top of the atmosphere.

But the term is large, and it is what makes the surface budget close against TOA:

| | TOA | SFC without snow | snow enthalpy | SFC (part2) | TOA − SFC |
|---|---:|---:|---:|---:|---:|
| AMIP 1872–79 | 0.667 | 1.441 | 0.925 | **0.516** | +0.152 |
| 09C 1370–79 | 0.918 | 1.897 | 0.944 | **0.953** | −0.035 |
| 080a 1370–79 | 1.682 | 2.588 | 0.823 | **1.765** | −0.083 |

The snow enthalpy is **0.82–0.94 W m⁻²** — larger than the entire remaining gap to
target. Without it the surface budget misses TOA by ~1 W m⁻². The 080a closure of
−0.083 reproduces the report's quoted ≈−0.09 exactly, confirming part2 is doing the
right thing and that the report's "surface tracks TOA to ≈0.1 W m⁻²" holds *only
because* the correction is applied.

**Convention warning.** The workplan's acceptance target is on the *surface*
(|F_sfc| ≲ 0.3 W m⁻²), not TOA. Surface-framed, the same runs read: atmosphere-only
floor **+0.52**, 09C **+0.95**, 080a **+1.77**. Always state which convention is in use.

*Open issue:* AMIP's closure (+0.152) is 2–4× worse than the coupled runs' (−0.035 /
−0.083). Small next to the +0.67 signal, so it does not threaten the attribution, but
it is unexplained. Likely candidates are the prescribed SST/sea-ice surface treatment
or an enthalpy term specific to the AMIP surface tile. Chase this before quoting the
AMIP floor to better than ~0.15 W m⁻².

---

## 3. The attribution result

**AMIP net TOA = +0.67 W m⁻²** (1872–79, after discarding two soil spin-up years;
per-year 0.13–0.92). `tsr` = 240.3 W m⁻² confirms the 3600 s divisor.

With SST prescribed, this cannot be ocean heat uptake — it is a pure atmospheric
radiative error. Against a coupled baseline of +1.44 (080a), **~45% of the coupled
imbalance is atmospheric**, and the best coupled branch (09C, +0.92) sits only
0.25 W m⁻² above the atmosphere-only floor. Of the 1.08 W m⁻² still separating 09C
from the HR goal (−0.16), ~0.83 is atmospheric.

*Caveat:* AMIP runs over *observed* SST/ice while the coupled model runs over its own
biased state, so +0.67 is "what this atmosphere does over a correct ocean", not a
strict lower bound on the coupled value. For a ~PI configuration where the real world
was near equilibrium the expected AMIP TOA is ≈0, so most of the +0.67 is model error.

---

## 4. Round 09 — new sea-ice thermodynamics

All 30 yr, last decade. 09A/09B are a270270's; 09C is ours.

| run | net TOA | global T2m | Siberia JJA | Siberia ann |
|---|---:|---:|---:|---:|
| 080a baseline, old sea ice | 1.442 | 13.606 | −1.41 | −1.36 |
| 09A baseline + newSeaIce | 1.286 | 13.148 | −2.19 | −3.00 |
| 09B 06T + newSeaIce | 1.171 | 12.794 | −2.50 | −3.11 |
| 08B 06V, old sea ice | 0.935 | 12.347 | −3.11 | −4.09 |
| **09C 06V + newSeaIce** | **0.919** | 12.492 | **−2.54** | **−3.08** |

On the clean old-vs-new 06V pair (08B → 09C) radiation barely moves (0.935 → 0.919)
but the boreal **warms by +1.0 K in the Siberian annual mean**, +0.57 K in JJA,
+0.15 K globally. 09C is the best radiation number in the campaign *and* a degree
warmer over Siberia than the previous best — it dominates 08B on both objectives.

The effect is **not** a simple offset: on the baseline (080a → 09A) the new ice
thermodynamics go the other way, cooling 0.46 K globally and 1.6 K in the Siberian
annual mean. The sign depends on the ocean-mixing branch. No mechanism established;
would need the surface energy budget over ice to pin down.

Coupling verified live via the OASIS `namcouple`: the 09 runs carry a sixth o2a field,
`sit_feom → A_Ice_thickness`, absent from all 08 runs.

---

## 5. Boreal seasonal bias — corrected AMIP build

model − CRUNCEP3 (1901–10), land only, cos-lat weighted.

| series | Boreal JJA/ann | Siberia JJA/ann | E. Siberia JJA/ann |
|---|---:|---:|---:|
| AMIP, **old** build (superseded) | −0.94 / −0.43 | −1.23 / −0.15 | −1.59 / +0.13 |
| **AMIP, corrected build** | −1.64 / −1.05 | **−2.16** / −0.89 | **−2.24** / −0.49 |
| 080a | −0.36 / −0.90 | −1.41 / −1.36 | −1.80 / −1.48 |
| 09A | −1.24 / −2.06 | −2.19 / −3.00 | −2.47 / −2.89 |
| 09B | −1.56 / −2.33 | −2.50 / −3.11 | −2.72 / −3.04 |
| 09C | −1.93 / −2.64 | −2.54 / −3.08 | −2.69 / −3.05 |

The corrected build is **substantially colder** than the old forcing run. The
atmosphere alone is now colder in summer than the coupled baseline and accounts for
essentially all of 09A's summer bias and most of 09B/09C's. **The boreal summer cold
bias is atmospheric.** The cold-season half is unchanged: AMIP's Siberian annual bias
is only −0.89 against −3.0 coupled, so Sep–Apr remains a coupled/vegetation problem.

**Data-quality defect:** the "Boreal 55-70N" box includes Greenland (JJA snow depth
320 mm w.e. in that box vs 5.5 mm for Siberia). The Siberia and E. Siberia boxes are
clean; the Boreal box is ice-sheet contaminated and should be masked or dropped.

---

## 6. Why the summer bias exists — JJA surface energy budget

Siberia JJA, land-masked, cos-lat weighted (AMIP 1872–79 vs CRUNCEP3 1901–10):

| | AMIP | CRUNCEP3 | Δ |
|---|---:|---:|---:|
| T2m | 9.74 °C | 11.92 °C | **−2.18 K** |
| SW net at surface | 155.1 | 188.3 | **−33.2 W m⁻²** |
| LW net at surface | −49.2 | −64.5 | +15.2 W m⁻² |
| net radiation | | | **−18 W m⁻²** |

Diagnostics from AMIP alone: SW↓ 188.6, surface albedo 0.178 (`fal` 0.172), LW↓ 315.6,
SH/LH −24.6/−56.8 (**Bowen 0.43**), cloud tcc 0.778 / lcc 0.505 / mcc 0.386,
skt−2t +0.43 K, `swvl1` 0.343, snow depth 5.5 mm w.e.

The LW term is compensation, not cause. **It is not albedo** — driving albedo to zero
recovers only 33 W m⁻², and a realistic 0.13–0.15 buys 5–9. So SW↓ itself is too low:
too much cloud, `lcc` 0.51 over Siberian land in high summer.

### The constraint problem

Burning off cloud returns SW to the surface but stops reflecting it at TOA. The Siberia
box is ~2.5% of the globe and JJA is 25% of the year, so +30 W m⁻² there is roughly
**+0.18 W m⁻² on the annual global TOA** — a ~27% worsening of the +0.67. The extra
OLR from a 2 K warmer surface returns only ~0.03. **A pure cloud fix does not meet a
fixed-TOA constraint.**

### The knob that does: surface flux partitioning — free at TOA

`susveg_mod.F90:264-267`, with the originals still commented out in place:

```fortran
!RVRSMIN(3)=500._JPRB    ! Evergreen Needleleaf Trees
!RVRSMIN(4)=500._JPRB    ! Deciduous Needleleaf Trees
 YDVEG%RVRSMIN(3)=250._JPRB    ! <- active
 YDVEG%RVRSMIN(4)=250._JPRB    ! <- active
```

Minimum stomatal resistance for boreal needleleaf — including larch — has been
**halved from the IFS default**. Half the resistance is double the conductance: more
transpiration, more latent heat, less sensible heat, colder screen temperature. That
matches the measured Bowen ratio of 0.43 (boreal conifer forest is typically 0.5–1.5).
Moving energy between SH and LH **changes no TOA radiation at all**; the only TOA
effect is the small OLR increase from a warmer surface, which *reduces* the imbalance.

### Proposed matrix — 5 yr each on AMIP, ~50 min per run

| # | change | where | expected |
|---|---|---|---|
| **A** | `RVRSMIN(3,4)` 250 → 500 (restore IFS default) | `susveg_mod.F90:266-267` | +1–2 K Siberia JJA, ΔTOA ≈ 0 |
| B | `RVLAI(3,4)` 5.0 → 4.0 | `susveg_mod.F90:133-134` | same sign, smaller |
| C | `RCLDIFF` ↑ (cloud erosion; NAMCLDP, namelist-settable) | runscript | +surface SW, but +TOA |
| D | A + C balanced so ΔTOA ≈ 0 | | joint target |

Run A first: the only one free at TOA, and it restores a documented default rather than
inventing a tuning.

**Caveats.** The Bowen-ratio claim is *inferred*, not measured — CRUNCEP3 carries no
SH/LH; FLUXNET boreal sites would settle it. Part of the low Bowen ratio is a
*consequence* of the SW deficit (less available energy shrinks both fluxes), so A alone
may not deliver the full 2 K and C may be needed regardless. CRUNCEP3 is a
reanalysis-derived forcing product, not observations. Period mismatch 1872–79 vs
1901–10 throughout.

---

## 7. Cost, for planning the tuning loop

AMIP: 10 yr in 1:38–2:45 on 14 nodes ≈ **2.3 node-hours per model year**.
Coupled: ~3.5 h per decade on 45 nodes ≈ **15.8 node-hours per model year**. ~7× cheaper.

Interannual scatter of the AMIP boreal JJA mean (measured, `iav.py`): σ = 0.25 K
(Boreal), 0.59 K (Siberia), 0.72 K (E. Siberia). With 8 usable years the standard error
is 0.09–0.26 K against biases of 1–2 K. Sampling is not the binding constraint —
**land spin-up is**: discard the first 1–2 years (`lresume: false`, deep soil
`stl4`/`swvl4`). For paired tuning experiments 5 yr (1 discard + 4 used) suffices,
since identical prescribed SST makes most interannual variance cancel in the difference.
