# Coupled 11E — bringing the coupled configuration level with AMIP

**For the colleague running the `Tuning_test_11*` series. Dated 2026-08-09.**

11D is the closest the coupled branch has ever been to the atmosphere-only configuration,
but it differs from it in two ways, and one of them matters.

## What 11D already has right

`ECE_TUNE_RVRSMIN` 1000/1000/225 (G4), `RCL_INPSEA` 0.2 with `RCL_INPPMIN` 70000 (D2b),
and the fitted snow depletion `ECE_SNOW_SCF=3`. That is most of the AMIP stack, and the
report was out of date in saying the coupled branch carried only G4.

## Change 1 — `SWEMIN` 30 → 15. This is the one that matters.

```yaml
                ECE_SNOW_SCF_SWEMIN: 15.0     # was 30.0
```

**Why.** `SWEMIN=30` is the P6 configuration, and P6 was *rejected* in AMIP over 44 years:

| AMIP run | DJF | MAM | JJA | SON | Nov `f_full` | DJF soil |
|---|---:|---:|---:|---:|---:|---:|
| P5 `SWEMIN=15` | **−0.060** | −0.140 | +1.149\* | +0.004 | 0.87 | **+0.22** |
| P6 `SWEMIN=30` | **−0.850\*** | −0.288 | +1.096\* | +0.487\* | **0.752** | **−1.00** |

\* clears that season's own 44-yr threshold (DJF ±0.588, MAM ±0.386, JJA ±0.242, SON ±0.431).

P6 lands on the "warms JJA but cools DJF significantly" list — the exact failure mode this
whole snow-scheme line has been trying to escape since round 15. The mechanism: the mass
floor binds in **November**, when the pack is still shallow, and November cover falls to
0.752 against the scheme-off 0.884. The summer gain is not worth it either — P6 buys
+0.060 K of JJA over the base, which is noise, against P5's +0.113.

**11D does not disprove this.** Its Siberian DJF soil is −0.04 K against 11A, but the
measured detection threshold on a 10-year coupled pair is **±0.75 K for soil and ±1.25 K
for DJF T2m** (from the interannual scatter of these very runs). P6's AMIP failure was
−1.00 K and −0.850 K. So the T2m half of the falsifier is *untestable* at this run length
and a null was the expected outcome either way. The diagnostic is not weak in general — it
sees 11B's −18.99 K without difficulty — it is just not able to resolve 1 K.

Take the AMIP answer. It says 15.

## Change 2 — add K1, the snow-free land albedo correction

Absent from 11D entirely. It is namelist-only, goes in the same `NAMSURFTUNE` block, and
is worth +0.089 K of global land T2m and the campaign's second-best global T2m RMSE
(1.524 against a 1.579 control).

```yaml
                # --- K1: vegetation albedo, types 1 / 10 / 11 only ---
                "ECE_TUNE_RVVEGALB(1,1)": 0.023500
                "ECE_TUNE_RVVEGALB(1,2)": 0.030240
                "ECE_TUNE_RVVEGALB(1,3)": 0.252561
                "ECE_TUNE_RVVEGALB(1,4)": 0.284663
                "ECE_TUNE_RVVEGALB(10,1)": 0.022748
                "ECE_TUNE_RVVEGALB(10,2)": 0.026933
                "ECE_TUNE_RVVEGALB(10,3)": 0.234026
                "ECE_TUNE_RVVEGALB(10,4)": 0.253862
                "ECE_TUNE_RVVEGALB(11,1)": 0.047943
                "ECE_TUNE_RVVEGALB(11,2)": 0.078395
                "ECE_TUNE_RVVEGALB(11,3)": 0.185338
                "ECE_TUNE_RVVEGALB(11,4)": 0.262827
                # --- K1: bare-soil background scale, all four components ---
                "ECE_TUNE_SOILALB(1)": 0.95
                "ECE_TUNE_SOILALB(2)": 0.95
                "ECE_TUNE_SOILALB(3)": 0.95
                "ECE_TUNE_SOILALB(4)": 0.95
```

Types 1 (crops), 10 (irrigated crops) and 11 (semidesert) only. **Tundra was deliberately
excluded** so that K stays orthogonal to the boreal levers and adds no absorbed shortwave
to the Arctic — 60–90N net TOA stays at −97.463 against a control −97.737 and CERES
−97.98. Do not add tundra entries: they compete for the same Arctic energy budget.

## Two traps, both paid for already

**`NAMSURFTUNE`, not `NAMECECFG`.** `NAMECECFG` is read twice from the same `fort.4` — by
`ECE_CONFIG` in `arpifs/ecearth.F90` and by `SURFECE_CONFIG` in `surf/surfece.F90` — and a
Fortran namelist read dies with *"invalid reference to variable"* on any name the reading
module does not declare.

**Build from a commit that has the `SWEMIN` floor.** Anything before `5952bed` hits the
`ABOR1 'Very snow cold temperature'` crash that killed P1 and 11C: the fitted power law
extrapolates below any depth the snow courses sampled, a 0.56 mm dusting gets full cover,
and the snow tile takes the whole grid-box flux into near-zero heat capacity. Also make
sure the tree carries `562df81`, which fixes `DCMAX` clipping the mass floor — without it
`SWEMIN` above 30 at low density is silently inert, which would make this whole change a
no-op at exactly the wrong end.

## What 11E would settle

Three things AMIP structurally cannot price, because it holds SST and sea ice fixed:
whether the boreal closure suffices for forest survival once the vegetation–albedo runaway
is allowed to close; whether G4's delay of melt-out from 24 to 28 May *harms* forest
establishment through later growing-season onset (a real risk, and AMIP is blind to it);
and the Arctic energy bill, which prescribed SST and ice absorb without responding.

If the run can be longer than 10 years, it is worth it — see the detection thresholds
above. Twenty years would roughly halve them and bring a 1 K DJF question into range.
