# OpenIFS 48r1 upstreaming: merge request titles and descriptions

The campaign branch `movcav-landice+co2-concdriven` split into seven stacked branches on `git.smhi.se/jan.streffing/oifs48r1`, to be merged into `ec-earth/vendor/openifs/oifs48r1` `main` in the order below. Each branch is cut from the one before it, so a branch's own change is its diff against its predecessor, and its diff against `main` shrinks to that once its predecessors have landed.

The "standalone" column records whether the branch's own commit also applies cleanly straight onto `main`. The first two do, so they can be reviewed and merged in either order. The rest overlap their predecessors textually, almost always in `surfece.F90`, and need the stack order.

Merge request numbers assume the six still to be opened take !92 through !97 in order, which holds only if nobody else opens one on this project first. !91 is confirmed. Check the numbers before pasting, and fix the stack sentence at the top of any description whose number came out different.

Copy the fenced block for each merge request straight into the GitLab description field. The fences are not part of the text: everything inside one is the description, already in GitLab-flavoured markdown.

| # | MR | branch | own diff | standalone |
|---|----|--------|----------|------------|
| 1 | !91 | `pr1/landice-ism-coupling` | 8 files, +198 / -35 | yes |
| 2 | !92 | `pr2/coupling-fixes-and-co2-spinup` | 5 files, +63 / -9 | yes |
| 3 | !93 | `pr3/ocean-skin-and-surface-fixes` | 10 files, +473 / -14 | no |
| 4 | !94 | `pr4/snow-cover-depletion` | 2 files, +502 / -7 | no |
| 5 | !95 | `pr5/dms-marine-ccn` | 7 files, +143 / -6 | no |
| 6 | !96 | `pr6/namelist-tuning-exposures` | 16 files, +432 / -11 | no |
| 7 | !97 | `pr7/lpjg-raupach-roughness` | 4 files, +194 / -14 | no |


---

## 1. `pr1/landice-ism-coupling`  (!91)

Commit `6eb10e1`, author Jorjo Bernales.

Already open as !91. Replace its title and description with the text below, because GitLab does not re-read either when a branch is force-pushed.

**Title**

```
Land-ice coupling: ice-sheet mask, ice surface physics and the ISM exchange
```

**Description**

```markdown
Part of an ordered stack of 7 merge requests, !91 through !97, splitting the AWI-ESM3 tuning campaign branch. This is the first of them, and it also applies cleanly to `main` on its own.

Adds an optional interactive land-ice surface to HTESSEL, together with the OASIS field exchange with an ice-sheet model (CISSEMBEL). Everything is gated off by default behind ECE_LANDICE and ECE_ISM_OROG, so an unconfigured run is bit-identical.

Over cells flagged as ice sheet, the soil column is replaced by ice thermodynamics and snow is coupled to an ice substrate. Roughness lengths are set for ice, as are the SW albedo and the LW emissivity, and the SEB linearisation is adjusted to match. The mask is either read from the forcing file or delivered by OASIS as PLIT, with a configurable threshold and fractional blending between ice and land physics. Multiple ISM regions are supported through a region loop in the coupling driver.

Squashed from 40 commits on the AWI-ESM3 tuning campaign branch. The originals
are listed at https://git.smhi.se/jan.streffing/oifs48r1/-/commits/movcav-landice%2Bco2-concdriven
and each can be read on its own at https://git.smhi.se/jan.streffing/oifs48r1/-/commit/<sha>:

- [`09f4de9`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/09f4de9) Define coupling stages/fields for OIFS-->ISMM (Violet Patterson)
- [`2382b71`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/2382b71) Add ISM region loop (Violet Patterson)
- [`e1d0708`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/e1d0708) tidy up (Violet Patterson)
- [`930bc84`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/930bc84) LandIce: Couple snow to ice substrate (Jorjo Bernales)
- [`0956bfb`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/0956bfb) LandIce: Replace soil by ice thermodynamics (Jorjo Bernales)
- [`bb99c16`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/bb99c16) LandIce: Impose ice roughness lengths (Jorjo Bernales)
- [`70a0be4`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/70a0be4) Tidy up (Violet Patterson)
- [`94447a4`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/94447a4) ISMM-->OIFS field exchange (Violet Patterson)
- [`d3bf27b`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/d3bf27b) LandIce: Modify SEB linearisation (Jorjo Bernales)
- [`be849dd`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/be849dd) LandIce: Adjust SW albedo and LW emissivity (Jorjo Bernales)
- [`a383624`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/a383624) LandIce: Implement evolving albedo over ice (Jorjo Bernales)
- [`41334b0`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/41334b0) LandIce: Add storage/API in SURFECE module (Jorjo Bernales)
- [`c3b2608`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/c3b2608) LandIce: Read/pack ice sheet mask from file (Jorjo Bernales)
- [`748c647`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/748c647) LandIce: Add active-block bridge for SURF routines (Jorjo Bernales)
- [`c7234e9`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/c7234e9) LandIce: Wire ice-sheet mask to surf routines (Jorjo Bernales)
- [`f535186`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/f535186) Add ECE_LANDICE gate and conditional mask wiring (Jorjo Bernales)
- [`d8ce225`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/d8ce225) Add ECE_LANDICE to ECEARTH module as well (Jorjo Bernales)
- [`8cbc4fc`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/8cbc4fc) LandIce: Configurable mask threshold and fractional physics blending (Jorjo Bernales)
- [`98fd0cf`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/98fd0cf) Framework for multiple ism regions (Violet Patterson)
- [`655097d`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/655097d) Fix mis-tag in DISGRID_RECV for PLIT read-in (Jorjo Bernales)
- [`4ef058a`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/4ef058a) Add region loop to Landice block (Violet Patterson)
- [`a0ec0ae`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/a0ec0ae) LandIce: Harden mask intake with error handling and diagnostics (Jorjo Bernales)
- [`06ff388`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/06ff388) LandIce: Fix thread-safety and consolidate surface physics constants (Jorjo Bernales)
- [`5383f36`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/5383f36) LandIce: Unified region loop for orography and mask reads (Jorjo Bernales)
- [`0e4032d`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/0e4032d) LandIce: Wire OASIS-delivered PLIT into SURFECE monthly update (Jorjo Bernales)
- [`21387d6`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/21387d6) LandIce: Drop duplicate ECE_LANDICE import in ece_updclie_cpl (Jorjo Bernales)
- [`393c21d`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/393c21d) LandIce: Read plit_oifs and usurf_oifs from forcing file (Jorjo Bernales)
- [`415091b`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/415091b) LandIce: Remove noisy unset-block-index warning from SURFECE (Jorjo Bernales)
- [`69315a3`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/69315a3) LandIce: Clean up development artifacts (Jorjo Bernales)
- [`e3f2efd`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/e3f2efd) LandIce: narrow USE NETCDF to ONLY used symbols in suorog (Jorjo Bernales)
- [`b9ba5df`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/b9ba5df) Restore original formatting in modified SURF routines (Jorjo Bernales)
- [`7b730d9`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/7b730d9) [OIFS] CISSEMBEL field exchange (Violet Patterson)
- [`7d043b2`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/7d043b2) [OIFS] allocate ISMRGNCIS (Violet Patterson)
- [`722fe38`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/722fe38) file merge edit corrections (Violet Patterson)
- [`a331739`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/a331739) Fix merge error (Violet Patterson)
- [`5779a3d`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/5779a3d) Create separate ECE_ISMM_SET_STATE subroutine, call in callpar.F90 (Violet Patterson)
- [`38c3ab0`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/38c3ab0) callpar: give the ECE_ISMM_SET_STATE call an explicit interface (Jan Streffing)
- [`83c720c`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/83c720c) Make the OASIS-driven land-ice mask update opt-in (Jan Streffing)
- [`460015f`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/460015f) suorog: gate the orography read on ECE_ISM_OROG (Jan Streffing)
- [`e13e341`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/e13e341) updtim LANDICE: describe what the hook does (Jan Streffing)

Co-authored-by: Violet Patterson <vip@dmi.dk>
Co-authored-by: Jan Streffing <jan.streffing@awi.de>
```

---

## 2. `pr2/coupling-fixes-and-co2-spinup`  (!92)

Commit `c5bc7ec`, author Jan Streffing.

**Title**

```
Coupling fixes and a concentration-driven CO2 spin-up mode
```

**Description**

```markdown
Part of an ordered stack of 7 merge requests, !91 through !97. It is branched from !91, but it also applies cleanly to `main` on its own, so it does not have to wait for !91 to land.

Three coupling changes, unrelated to each other:

- initialise OASIS_OUT fields to 0 rather than -HUGE, so the first send does not carry the sentinel into the receiving model
- receive FESOM sea-ice thickness for the coupled-slab ice surface
- add ECE_CPL_FESOM_RECOM_CONCDRIVEN, a concentration-driven CO2 spin-up mode for FESOM-RECOM

Squashed from 3 commits on the AWI-ESM3 tuning campaign branch. The originals
are listed at https://git.smhi.se/jan.streffing/oifs48r1/-/commits/movcav-landice%2Bco2-concdriven
and each can be read on its own at https://git.smhi.se/jan.streffing/oifs48r1/-/commit/<sha>:

- [`b7eacd0`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/b7eacd0) cplng2: init OASIS_OUT fields to 0, not -HUGE (fixes first-send sentinel) (Jan Streffing)
- [`058c473`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/058c473) fesom-fesim: coupled-slab ice surface -- receive FESOM ice thickness (Jan Streffing)
- [`d3e3817`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/d3e3817) Add ECE_CPL_FESOM_RECOM_CONCDRIVEN spinup mode for CO2 (Jan Streffing)
```

---

## 3. `pr3/ocean-skin-and-surface-fixes`  (!93)

Commit `da5ce09`, author Jan Streffing.

**Title**

```
Ocean skin: only over water, and fix warm-start surface field seeding
```

**Description**

```markdown
Part of an ordered stack of 7 merge requests, !91 through !97. It is branched from !92 and overlaps its predecessors textually, so !91 through !92 need to land first. Until they do, the diff shown here against `main` includes their changes as well.

The mesh-increment friction velocity could overflow where voskin was evaluated over cells with no water. It is now computed only where there is ocean, and the diagnostic probes added to find that are removed again.

Separately, two warm-start defects. reresf did not seed surface defaults before re-reading ICMGG, and surface_fields_mix SETDEFAULT skipped fields with NREQIN==0. Both left surface fields undefined on a warm start. The 4-layer ice temperature now also uses the coupled ocean sea-ice thickness.

Squashed from 8 commits on the AWI-ESM3 tuning campaign branch. The originals
are listed at https://git.smhi.se/jan.streffing/oifs48r1/-/commits/movcav-landice%2Bco2-concdriven
and each can be read on its own at https://git.smhi.se/jan.streffing/oifs48r1/-/commit/<sha>:

- [`a284bd1`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/a284bd1) OIFS moving-cavity/moving-lsm coupling fixes (Jan Streffing)
- [`491fb76`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/491fb76) surf: use coupled ocean sea-ice thickness in the 4-layer ice temperature (Jan Streffing)
- [`f3ccacb`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/f3ccacb) fix surfece switch (Jan Streffing)
- [`45a7bea`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/45a7bea) voskin: diagnostic probes for the mesh-increment friction-velocity overflow (Jan Streffing)
- [`3a918c8`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/3a918c8) reresf: seed surface defaults before re-reading ICMGG on a warm start (Xiaojie Hao)
- [`ab9f53b`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/ab9f53b) surface_fields_mix: SETDEFAULT also seeds NREQIN==0 fields (Xiaojie Hao)
- [`04d1481`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/04d1481) voskin: only compute the ocean skin where there is water (Xiaojie Hao)
- [`f00a397`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/f00a397) strip the diagnostics the voskin root cause made redundant (Xiaojie Hao)

Co-authored-by: Xiaojie Hao <xiaojie.hao@awi.de>
```

---

## 4. `pr4/snow-cover-depletion`  (!94)

Commit `9ec874b`, author Jan Streffing.

**Title**

```
Optional snow-cover depletion formulations (ECE_SNOW_SCF), off by default
```

**Description**

```markdown
Part of an ordered stack of 7 merge requests, !91 through !97. It is branched from !93 and overlaps its predecessors textually, so !91 through !93 need to land first. Until they do, the diff shown here against `main` includes their changes as well.

Adds alternative snow-cover fraction parameterisations, selected by ECE_SNOW_SCF and all off by default:

- 1: unchanged operational behaviour
- 2: scale-aware depletion via SDOR
- 3: depletion fitted to observations, resolution-aware

Mode 3 floors the depletion depth at a minimum snow mass (ECE_SNOW_SCF_SWEMIN) rather than at a depth, which fixes a crash where DCMAX could otherwise disable the mass floor entirely. ECE_SNOW_SCF_Z0 defaults to the calibrated 0.016 m. Permanent-snow resets are extended to thin-snow glacier cells in AMIP configurations.

Squashed from 7 commits on the AWI-ESM3 tuning campaign branch. The originals
are listed at https://git.smhi.se/jan.streffing/oifs48r1/-/commits/movcav-landice%2Bco2-concdriven
and each can be read on its own at https://git.smhi.se/jan.streffing/oifs48r1/-/commit/<sha>:

- [`0c9e575`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/0c9e575) Extend permanent-snow resets to thin-snow glacier cells in AMIP (Jorjo Bernales)
- [`28b5542`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/28b5542) Add optional Niu & Yang (2007) snow-cover depletion, off by default (Jan Streffing)
- [`9f12d79`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/9f12d79) Set ECE_SNOW_SCF_Z0 default to the calibrated 0.016 m (Jan Streffing)
- [`6e817ea`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/6e817ea) Add ECE_SNOW_SCF=2: scale-aware depletion via SDOR (Jan Streffing)
- [`c16220a`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/c16220a) Add ECE_SNOW_SCF=3: snow-cover depletion fitted to observations (Jan Streffing)
- [`5952bed`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/5952bed) Fix ECE_SNOW_SCF=3 crash: floor d_c at a minimum snow MASS (ECE_SNOW_SCF_SWEMIN) (Jan Streffing)
- [`562df81`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/562df81) Make ECE_SNOW_SCF=3 resolution-aware, and stop DCMAX from disabling the mass floor (Jan Streffing)

Co-authored-by: Jorjo Bernales <bernales.work@gmail.com>
```

---

## 5. `pr5/dms-marine-ccn`  (!95)

Commit `c858bdb`, author Jan Streffing.

**Title**

```
Marine biogenic CCN from DMS, with the coefficient in the namelist
```

**Description**

```markdown
Part of an ordered stack of 7 merge requests, !91 through !97. It is branched from !94 and overlaps its predecessors textually, so !91 through !94 need to land first. Until they do, the diff shown here against `main` includes their changes as well.

Wires a DMS surface concentration field through the ICMCL ingest chain into PCCNO, so that marine biogenic aerosol can modulate the cloud droplet number over ocean instead of the fixed climatological value. The conversion coefficient is exposed as a namelist parameter, and with no DMS field present the behaviour is unchanged.

Squashed from 1 commit on the AWI-ESM3 tuning campaign branch. The originals
are listed at https://git.smhi.se/jan.streffing/oifs48r1/-/commits/movcav-landice%2Bco2-concdriven
and each can be read on its own at https://git.smhi.se/jan.streffing/oifs48r1/-/commit/<sha>:

- [`8707cd1`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/8707cd1) Wire marine biogenic CCN from DMS into PCCNO, with the coefficient as a namelist knob (Jan Streffing)
```

---

## 6. `pr6/namelist-tuning-exposures`  (!96)

Commit `bf38056`, author Jan Streffing.

**Title**

```
Expose HTESSEL, cloud and stable-BL tuning constants to namelists
```

**Description**

```markdown
Part of an ordered stack of 7 merge requests, !91 through !97. It is branched from !95 and overlaps its predecessors textually, so !91 through !95 need to land first. Until they do, the diff shown here against `main` includes their changes as well.

Makes tuning constants that were hard-coded PARAMETERs settable from the namelist, with every default left at its as-released value, so this changes no results by itself. It exists so that a tuning campaign does not need a recompile per parameter.

- NAMSURFTUNE: the HTESSEL tuning tables, the snow-tile skin conductivities and RVLAMSKS
- NAMCLDP: RCLCRIT_SEA / RCLCRIT_LAND and two parameters for continental INP, plus RCAPDCYCL mode 4
- NAMVDFS: the stable-boundary-layer mixing constants (RSBLB, RSBLPBLF, RSBLLMAX and the length-scale floor)

Squashed from 7 commits on the AWI-ESM3 tuning campaign branch. The originals
are listed at https://git.smhi.se/jan.streffing/oifs48r1/-/commits/movcav-landice%2Bco2-concdriven
and each can be read on its own at https://git.smhi.se/jan.streffing/oifs48r1/-/commit/<sha>:

- [`5bea52e`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/5bea52e) Round 10 tuning-campaign annotations (values at as-released defaults) (Jan Streffing)
- [`2630bd1`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/2630bd1) Round 11: RCAPDCYCL mode 4, and two NAMCLDP parameters for continental INP (Jan Streffing)
- [`1004cba`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/1004cba) Make the HTESSEL tuning tables namelist-settable via NAMSURFTUNE (Jan Streffing)
- [`509c7b0`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/509c7b0) Expose the snow-tile skin conductivities to NAMSURFTUNE (Jan Streffing)
- [`068d523`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/068d523) Expose RVLAMSKS to NAMSURFTUNE (Jan Streffing)
- [`75a4e45`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/75a4e45) Complete the DMS ICMCL ingest chain, and expose RCLCRIT_SEA/LAND (Jan Streffing)
- [`5893ee0`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/5893ee0) Expose the stable-BL mixing constants to NAMVDFS (Jan Streffing)
```

---

## 7. `pr7/lpjg-raupach-roughness`  (!97)

Commit `9abf1fa`, author Jan Streffing.

**Title**

```
Receive the Raupach canopy roughness length from LPJ-GUESS
```

**Description**

```markdown
Part of an ordered stack of 7 merge requests, !91 through !97. It is branched from !96 and overlaps its predecessors textually, so !91 through !96 need to land first. Until they do, the diff shown here against `main` includes their changes as well.

When coupled to LPJ-GUESS, take the momentum roughness length for vegetated tiles from the dynamic vegetation's Raupach canopy formulation instead of the static per-type table, so that z0m follows the simulated canopy structure. Falls back to the table when the field is absent.

Squashed from 1 commit on the AWI-ESM3 tuning campaign branch. The originals
are listed at https://git.smhi.se/jan.streffing/oifs48r1/-/commits/movcav-landice%2Bco2-concdriven
and each can be read on its own at https://git.smhi.se/jan.streffing/oifs48r1/-/commit/<sha>:

- [`869be00`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/869be00) Receive the Raupach canopy roughness from LPJ-GUESS (Jan Streffing)
```
