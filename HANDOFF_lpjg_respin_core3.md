# Handoff: rebuild the LPJ-GUESS forcing and spin-up for the new CORE3 mesh

Written 2026-08-26, updated the same evening with the state of play. Everything below was
established today and is worth not re-deriving.

## Why this is needed

The 40-year mesh test `11X` (new CORE3, otherwise identical to `11W`) crash-loops. The
failure is **not** the mesh and not the atmosphere. LPJ-GUESS ranks exit 99 on two
independent checks that are new in tag **4.1.11** (`8832661`, Laszlo Hajdu, 2026-08-24):

```
Fixed peat invariant failed in after coupled restart load at lon 54.782608, lat 72.466919,
year 1350: map=0.0086111099999999999 physical=0 LC=0 ST=0; totals physical=1 LC=1 ST=1.
State was not repaired.

LPJG: configured restart could not be loaded for lon 61.578949, lat 76.207008; refusing
to cold-start this grid cell
```

1. `validate_peatland_invariant` (`framework/externalinput.cpp:4369-4457`, called from
   `framework.cpp:4481` after every restart load) compares the peat map against the stands in
   the state. It **validates only; there is no repair and no migration path** for a state that
   predates fixed-peat handling. It is gated solely by `run_landcover && run[PEATLAND]`, both
   on. "State was not repaired" is the fixed suffix of the message, not a failed attempt.
2. `framework.cpp:4448`: with `run[PEATLAND]` on, a gridlist cell with no entry in the state
   is fatal. The cold-start path of commit `1481b4f` only runs with peatland off.

**So the blocker is the STATE, not the code, and it is independent of the mesh.** The
state in use, `/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_2000Y_TCO95_CORE3_CRUNCEPandCERES_daily_variability/run_19000101-38991231/work/lpjg_state_3901`,
was made 2026-07-29 by a binary built 2026-07-07 (no fixed-peat handling), forced by
CRUNCEP + CERES, on an **11538-cell** gridlist derived from `L096.msk`. The old binary
`guess.raupach_junebase_gated` reads it fine, which is how `11W` ran 50 years to 1400.

The new mesh compounds it: the coastline-following gridlist (`TCO95-land.msk` of the new
`masks.nc`) has **11313** cells, so LPJ-GUESS meets cells the spin-up never simulated (the
second message above). Check 2 makes that fatal on its own.

Relaxing either check is the wrong move: this campaign has already decided once that
widening a tolerance does not restore comparability. `state_remap 1` (nearest-neighbour
state remapping, `ecearth.ins`) would only address check 2.

The "cheap discriminator" of running 4.1.11 on the old mesh with the old state was **not
run**: the code path settles it. Check 1 reads the same cell of the same state against the
same peat map on either mesh, so it fires identically on `core3_beta`.

## The job, and where it stands

**1. New 10-year AMIP run** on the adopted atmospheric stack, writing the daily LPJ-GUESS
forcing fields. **DONE 2026-08-26: `amip_LF2_lpjgforce`, job 27268698, COMPLETED in 1h35, ten
`outdata/oifs/atm_1d_18??-18??.nc` of ~560 MB.** Runscript
`~/esm_tools/runscripts/oifsamip/oifsamip-cy48-levante-TCO95L91_LF2_lpjgforcing_sb2nos4.yaml`.
Staged `fort.4` verified: `RSBLB=2`, `RCL_INPPMIN=70000`, `NCMIPFIXYR=1850`, `GGAUSSB=-0.5`,
`ENTSTPC3=1`, `RSNOWLIN2=0.04`; `field_def.xml` divides fluxes by 3600; staged
`libarpifs.SP.so` md5 `b3e563e3…` = the AMIP tree's build library. ~2 h of compute.
Year 1870 closed after ~11 min with sane magnitudes in `work/atm_1d_1870-1870.nc`:
global-mean `rsns` 164 W/m², `rlns` −60, `pr` 3.2e-5 kg m⁻² s⁻¹, so the 3600 s
de-accumulation is right (the 21600 bug would read ~27 W/m²).

Forcing file **DONE**: `scripts/sbatch/make_lpjg_forcing_LF2.sbatch` (job 27272265) wrote
`input/lpj-guess/oifs_forcing/AMIP_LF2_1d_1870-1879_TCO95_PI_sb2nos4.nc` (3.2 GB, 3652 days,
tas/tasmin/tasmax/rsns/rlns/pr/sfcWind/hurs; global means tas 287.3 K, rsns 164.4 W/m²,
pr 3.22e-5, hurs 74.7 % against the fluxfix file's 287.0 / 166.5 / 3.21e-5 / 75.2). (The generic one is `make_lpjg_forcing.sbatch`: edit `IN` to
`/work/bb1469/a270092/runtime/oifsamip-cy48/amip_LF2_lpjgforce/outdata/oifs`, `OUT` to
`/work/ab0246/a270092/input/lpj-guess/oifs_forcing/AMIP_LF2_1d_1870-1879_TCO95_PI_sb2nos4.nc`).
`create_lpjg_forcing.py` copies `rsns`/`rlns`/`pr` unchanged when XIOS already emitted them
(it only divides by 21600 for raw `ssr`/`str`/`tp`, which this file_def does not write), so
no post-hoc fluxfix is needed. **Check `rsns` over land is O(100) W/m²** before use.

**2. New 2000-year LPJ-GUESS spin-up** under 4.1.11 on the new mask. Runscript written:
`~/esm_tools/runscripts/lpj_guess/lpj_guess_levante_spinup_TCO95_CORE3_PIforcing_2000y_4111.yaml`.
It points at the forcing file above (now present). Submit only after the second smoke
test passes. Expect ~5 h at 1280 tasks (v3 took 4h54m). One job, cold start; chunked/lresume OOMs.

**Smoke test first**: `lpj_guess_levante_spinup_TCO95_CORE3_smoke_4111.yaml`, expid
`CORE3_smoke_4111`, 3 years on the existing fluxfix forcing, same mask/peat/LUH3 as
production. The 4.1.11 binary had never run in `-islpjgspinup` mode before this.

**First attempt (2026-08-26, job 27268932, kept as `CORE3_smoke_4111.failed_luh3_soilinit_20260826`)
looked healthy and was not.** 1276 of 1280 state files, output for every rank, gridlist
11313 cells as expected — but `scripts/analysis/lpjg_rank_outcomes.py` (parses every
`run<r+1>/guess<r>.log`) shows only **4563 of 10157 simulated cells finished**. In spin-up
mode each rank works through its cells sequentially and writes the state per cell, so a
rank that dies on its fourth cell leaves a plausible state file. Two independent killers:

1. **912 ranks**: `Cannot initialize soil for newly created land-cover stand type 0`
   (stand type 0 = `Natural`). 4.1.11 (`2de5d32`) added an early `soil_temp_multilayer`
   call in `Stand::init_stand_lu` so a stand created by land-cover change has its soil
   geometry before water is transferred in. On a gridcell's first simulated day that call
   runs inside `getgridcell`, before the cell's forcing is read and before
   `ifs_soil_parameters` (`framework.cpp:4473` vs `:4484`), so the solver returns false
   (`DIAG T_ERROR … T[7]=-4e21 … year=1900 day=0`) and `fail()` exits 99. Peat and cropland
   presence do not discriminate. Fix: commit `917c4a4` on the branch turns that `fail()`
   into a warning; `firstTempCalc` is still set, so the daily call on the same day
   initialises the soil as it always did before 4.1.11. Tell Laszlo.
2. **42 ranks**: `LUH3: all 12 required state fields are missing … refusing to substitute
   zero/NATURAL`. The coastline gridlist (`TCO95-land.msk`) contains cells LUH3 never
   covered: the Caspian and Aral, the Great Lakes, Victoria/Tanganyika/Malawi, a few Arctic
   islands and the Amazon mouth — 41 of them are water in the AMIP `slt`, i.e. lakes that
   ocp-tool now gives a soil type, so `ifs_soilcd > 0` no longer skips them. 4.1.11 requires
   an explicit audited list: `input/lpj-guess/land_use/TCO95/luh3_natural_override_TCO95_CORE3.txt`
   (47 cells = simulated cells with all 12 states at fill; Antarctica is skipped by
   `lat > -60` and is not in it), passed as `lpj_guess.file_luh3_natural_override`, now
   templated into both `guess.ins*.j2`. **11X needs this file too** — its logs already
   carry the same message.

Also seen: esm_tools' tidy of 1280 rank directories took >33 min for a 3-year run and hit
the 40-min limit; v3's "4h54m" was ~1 h of model and ~4 h of tidy. The state in `work/`
is complete before tidy starts; give the smoke 1h30 and production 8 h.

**Second attempt PASSED** (job 27272874, expid `CORE3_smoke_4111`, setup's own `bin/guess`
md5 `19169b1a…` = branch tip `917c4a4`, override file on): `lpjg_rank_outcomes.py` gives
1280/1280 ranks `ok`, all 10157 simulated cells `ok` (the 1156 "unreached" are Antarctica,
skipped by `lat > -60`), 10702 deferred-init warnings where the deaths were, no
"Cannot initialize"/"LUH3: all 12"/"Fixed peat invariant"/"refusing to cold-start". The
Caspian cell (49.09E, 42.55N) simulated as 100 % natural. Model step 2.3 min, job 1h04
(the rest is tidy), then esm_tools' 1-node `Combine_LPJG` subjob on the pp partition.

**Production spin-up**: expid `CORE3_2000y_LF2_4111`, **job 27273220**, submitted 2026-08-26 18:37; staged
`guess` md5 `19169b1a…`, `nyear_spinup 2000 / freenyears 100`, override on, forcing = the LF2 file (md5-checked). Expect ~1 h
of model plus hours of tidy inside the 8 h limit. Its product is
`run_19000101-38991231/work/lpjg_state_3901`. **Before using it**, run
`scripts/analysis/lpjg_rank_outcomes.py <that work dir> 1280` and require 1280/1280 ranks
`ok`; the state-file count proves nothing (see the first smoke attempt).

**Production PASSED** (2026-08-26 23:07): model step 33 min (22:34-23:07), `lpjg_rank_outcomes.py`
1280/1280 ranks `ok`, all 10157 cells, `work/lpjg_state_3901` = 1280 files + `meta.bin`,
3.4 GB, no fatal messages. Final-year (3899) vegetation carbon matches v3: mean `cmass`
5.19 vs 5.22 kgC/m², 67 % vs 64 % of cells above 1 kgC/m², so the biosphere is alive.
Table in `plots/lpjg_2000y_LF2_4111_cell_outcomes.csv`. The SLURM job was cancelled at
23:40 during tidy (1h06 elapsed; not by this session), so `restart/lpj_guess` and
`outdata/lpj_guess` are empty and everything stays in
`run_19000101-38991231/work/` — state (1281 files, 3.4 GB, intact) and per-rank output.
11X reads the state from that work directory, and its own staged copy is complete.

**11X is re-pinned and ready** (`awiesm3-develop-levante-TCO95L91-CORE3_11X_core3new_1850.yaml`):
`executable: guess.4111_landgrid_raupach_softinit` (coupled tree rebuilt from branch tip
`917c4a4` with `scripts/model/comp_lpjg_worktree.sh`, md5 `dc32c2fcfdbca91704871dbcff826bd3`;
the old `guess.4111_landgrid_raupach` is kept), `file_luh3_natural_override` set, and
`ini_restart_dir` pointing at the new spin-up. Resubmitted 2026-08-26 23:20 (job 27279252)
after the old directory was moved to `11X.failed_4111_oldstate_20260826`.

**That leg died at OpenIFS step 1689 (11 March 1350), not in LPJ-GUESS**: the state loaded,
none of the four LPJ-GUESS signatures appeared. `ABOR1: Very snow cold temperature`
(`srfsn_webals_mod.F90:452`, the `PTSN < 100 K` guard) on a 0.1 mm pack with snow fraction
0.008 — the same guard as the P1 crash in the report (§r20crash), but with the SWEMIN floor
already active (`fort.4` has `ECE_SNOW_SCF_SWEMIN=15`). Not an OIFS branch switch: the
checkout from `campaign` to `movcav-landice+co2-concdriven` was 08-23 09:44 and the
library 11W ran on was built 08-23 21:47 from the same commit `3c35275` as 11X's; the
tree is clean. Kept as `11X.failed_snowcold_nolakes_20260826`.

**The input defect behind it: the new-mesh ICMGG had no lakes.** `ICMGGab45INIT_CORE3_v2`
(ocp-tool regeneration of 2026-08-25, `logs/regen_core3_20260825.log`) has 0 land cells
with `cl > 0.5`; both old-mesh files have 70 (Caspian 27, Great Lakes 17, Aral, Victoria,
Tanganyika, Malawi, …). ocp-tool's `LakeConfig.restore_flipped_lakes` defaults to False on
the static path (`config.py:324`; only `dynamic_regen.py:67` forces it on), so every cell
flipped ocean→land inherited a dry neighbour's `cl` (0.0–0.45). This is also why the
LPJ-GUESS gridlist met those cells with a soil type. Repaired in place with
`ocp-tool/tools/repair_lake_cover.py` (pristine `ICMGGab45INIT`, published `_v2`) →
`ICMGGab45INIT_CORE3_v3`: cl/dl/FLake restored at the 102 flipped cells, `lsm` and `slt`
bit-identical (masks, gridlist, LPJ-GUESS soil file unchanged), read-back 68 land cells
with `cl ≥ 0.5` plus 2 Sea-of-Azov cells that are now FESOM ocean but keep `cl` (harmless
inconsistency, noted, not fixed). 11X runscript: `general.icmgg_suffix: "_v3"`.

**Resubmitted 2026-08-26 23:50 on `_v3`: job 27279956** (staged `ICMGGawi3INIT` verified: 11313 land, 68 lake cells). **Passed step 1689 and finished
model year 1350 in ~20 min** with LPJ-GUESS output flowing (`LPJ-GUESS_monthlyoutput.txt` 5.3 MB) and
the 1350 daily atmosphere files written, so the lake repair removed the abort; leg 1 (10 yr) ≈ 3.5 h. Watch leg 1 past step 1689 and to the point
where LPJ-GUESS output appears (`outdata/lpj_guess`), not the job state: a crashed leg still
reports COMPLETED and the date file advances. If the snow abort recurs on `_v3`, the lake
loss was not its cause and the point has to be located (only `JL` on OIFS task 1986 is
printed; the abort block would need lat/lon passed down from `surfexcdriver`).

## 2026-08-27: 11X on `_v3` died in FESOM, and the corrected-CORE3 coupled ocean is losing heat everywhere

Job 27279956 got past the snow abort, finished year 1350 with LPJ-GUESS coupling verified,
and was killed by esm_tools at OpenIFS step 9802 (≈12 Feb 1351) on a real
`MPI_ABORT` from FESOM rank 932: `found temperature becomes NaN or <-5.0, >60`, node
169679 (133.03E, 9.86S, Timor Sea, 100 m column, level 11 at −5.004 °C). Blow-up file
`run_13500101-13591231/work/fesom.1351.oce.blowup.nc`. Not a CFL event (cfl_z 0.09 there).

**It is not one node.** At the blow-up 5300 tropical nodes have SST < 10 °C and 735 carry
sea ice (e.g. 45W/23S off Santos at −2 °C, 30–80 % ice). They are the shallow shelves
(median depth 78 m vs 2884 m for all tropical nodes): Indonesia/Australia 1304, Brazil 476,
W Indian Ocean 470, Atlantic Africa 279. `core3y1` — the "sane" one-year zstar run of
08-26 — already had this: tropical SST<10 °C nodes 0 → 150 → … → 1513 over its year, and
its `fh` (heat flux, positive = ocean loss) averaged over the tropics is **+373 W/m² in
month 1** against **−8** in 11W's month 1 (11W: 40–60N +150, 60–90S −55, i.e. normal;
core3y1: 40–60N −1, 60–90S −34). `plots/` has nothing for this yet; numbers are from
`scripts` run inline, see the session transcript.

**What is right** (checked with the OASIS EXPOUT dumps left in
`runtime/awiesm3-develop/core3y1_expout_104302/run_*/work/`, mapped to mesh nodes by
matching FESOM's written `feom.lon/lat` to `nod2d.out`, which equals the rank-concatenated
`my_list` order): OpenIFS sends tropical Qns −188 / Qs +260 W/m²; FESOM receives −169 /
+236–333 at the right nodes (sent SST vs model SST corr 0.98); OASIS grids/masks/rmp are all
220509-node and normalised (weight sums 1.000, links median 227 km); A096 mask = OIFS lsm;
only 12 % of the cold nodes sit under atmosphere land. Received net heat in the tropics is
**+22.7 W/m² (gain)** while FESOM applies **+373 (loss)**, uncorrelated (r = −0.13).
So the damage happens inside FESOM between `oasis_get` and `heat_flux`.

**What is excluded**: an OIFS branch switch (11W and 11X both from `3c35275`);
`force_flux_consv` (returns at once under `__oifs`); the merged `__recom && __usetp` code
(`RECOM_COUPLED=OFF`); the tracer boundary condition (`bc_surface` unchanged for T);
partition numbering (APPLE partition is self-consistent). The 08-26 FESOM merge of main
`f4150b1d` (`4d7da170`) touched 41 files, 2769 lines, including `gen_surface_forcing.F90`
(860), `oce_ale_tracer.F90`, `fesom_module.F90`, `oce_setup_step.F90`, `io_restart.F90`;
11W ran the pre-merge `libfesom.so` (08-23, md5 `a1293118`), `core3y1`/11X the post-merge
one (md5 `39281712`).

**Decisive A/B submitted 2026-08-27 ~07:50**: `11Wab_newbin_oldmesh_1m` = 11W's exact
configuration (old mesh `core3_beta`, old state, old LPJ-GUESS binary) on today's binaries,
one month. Read `outdata/fesom/fh.fesom.1350.nc`, tropical mean of month 1: ≈ −8 ⇒ the
new-mesh inputs are at fault; ≈ +370 ⇒ the FESOM merge broke the coupling and the fix is
a bisection of `4d7da170` (or rebuilding from `backup/awiesm3-implicit-ice-surftemp_20260825`
plus #875 alone).

## The atmospheric stack for the forcing run (do not improvise)

Base **SB2** (LX4 + RSBLB), the adopted stack the coupled 11-series carries and what gave
the best campaign CMPI (11V, 0.7385). **Not SC3**: it is merely the newest AMIP file, and
`RDEPLIQREFRATE 0.2` is in the discarded table (SC3 breaks the Arctic, −1.85 on a band
already −18 over-reflective). A rejected lever in the forcing is inherited by every
downstream spin-up.

Kept: `RCL_INPSEA 0.2` · `RVICE 0.16` · `RSNOWLIN2 0.04` · `RSBLB 2.0` ·
`ECE_TUNE_RVRSMIN(3,4) 1000`, `(9) 225` · `ECE_SNOW_SCF 3` + its nine constants ·
`ECE_TUNE_RVVEGALB` 12 values · `SOILALB 0.95` · `ENTSTPC3 1` · `GGAUSSB −0.5` ·
`ECE_CLIMR_DMS true` with `ECE_DMS_CCN_SENS 0.0`. `LRDALB False` and the melt-pond settings
are coupled-only (LPJG albedo, FESOM) and do not apply to an AMIP run.

**S4 is retired, `RCL_INPPMIN = 70000`.** The review finished 2026-08-24: coupled, removing
S4 improves Siberian JJA by +0.590 K resolved, recovers the sub-Antarctic overcorrection
(+3.18/+3.45 W/m², replicated at both forcings), leaves the polar band untouched, and gives
the best CMPI. Every pre-existing AMIP runscript still carries 50000; none is a valid base.

**OIFS trees are separate.** AMIP builds from `model_codes/oifsamip-cy48`, coupled from
`model_codes/awiesm3-develop/oifs-48r1`. `&NAMVDFS` came from
`scripts/patches/expose_sbl_mixing.py` and had to be applied to both; verified today that the
AMIP tree carries `suvdfs.F90`/`vdfexcu.F90` with `RSBLB` and that the staged library is that
build. An absent group silently keeps 5.0.

## Facts you should not have to rediscover

**Mesh discipline.** `core3` = 220509 nodes, `core3_beta` = 211567. The shipped mesh was
replaced in place on 2026-08-26 and the old one renamed; the fix concerned cavity depths.
`restart/CORE3` no longer exists (`restart/core3_beta` holds the old `fesom.2000.*`;
`restart/core3_linfs` and `core3_zstar` the new-mesh 1959 restarts). Analyse pre-2026-08-26
runs with `core3_beta`, later runs with `core3`, never mix: numpy will area-weight a
211567-node field against the first 211567 entries of a 220509-node array and return a
plausible wrong number. Use `scripts/analysis/meshguard.py`, which raises instead. The new
`core3` ships **without** `fesom.mesh.diag.nc`; reval needs that file.

**Which inputs are new-mesh.** `input/oasis/cy48r1/TCO95-CORE3/masks.nc` md5 `fb15dc2d`
(old: `TCO95-core3_beta/`); `input/lpj-guess/slt/slt_TCO95_CORE3.nc` md5 `c00037f1` (=
ocp-tool `slt_TCO95_CORE3_v2.nc`, 2026-08-25); `peat/TCO95_peat_frac.txt` regenerated
2026-08-23 at full precision (`.orig2dp` is the old 2-decimal file). The lpjg-spinup setup
selects the first two through `fesom_resolution: CORE3`; the gridlist is generated from
`TCO95-land.msk` at run time, there is no gridlist file to install.

**Binary.** `model_codes/awiesm3-develop/bin/guess.4111_landgrid_raupach`, md5
`4b46d3564ab169ee52ccdbabfc5e485b`, from branch `feat/lpjg-land-grid-from-coastline-4.1.11`
= tag 4.1.11 plus the four coastline-branch commits (the `-land` gridlist, cold-start of a
cell with no state, the two Raupach `GUE_Z0HV` commits). The LUH3 tolerance workarounds and
the peat import/revert pair were dropped because 4.1.11 fixes the cause. **The branch is
local only, unpushed.** Build with `esm_master recomp-lpj_guess-4.1.2`, never cmake
directly, and verify the md5 changed. OASIS is linked statically, so the binary runs from the
`lpjg-spinup` setup; the spin-up runscripts point `bin_sources` at the coupled tree's `bin/`.

**The `lpjg-spinup-develop` tree was stale.** `model_codes/lpjg-spinup-develop/lpj_guess` was
the July `lpj_guess_awiesm3` checkout (`1cea929`), and its `bin/guess` (2026-07-05, md5
`89de1be8…`) is what any spin-up runscript without an explicit `executable`/`bin_sources`
stages. On 2026-08-26 evening the 4.1.11 coastline branch was fetched into it from the
coupled tree (remote `coupled`), checked out, and rebuilt with
`esm_master recomp-lpjg-spinup-develop`. The first attempt died in OASIS: that tree's OASIS
(upstream `main`, `a1ecb81`) hands `$ENV{OASIS_FFLAGS}` straight to `target_compile_options`,
and the esm_tools gnu block (`oasis3mct.env.yaml`, levante/gnu11_ompi4) exports five flags
as one string, which CMake never word-splits, so gfortran saw one unknown option. Patched the
four `lib/*/CMakeLists.txt` with `separate_arguments` (uncommitted, `git diff` in
`lpjg-spinup-develop/oasis`), wiped `oasis/build`, reran; log
`model_codes/comp_lpjgspinup_4111_20260826b.log`. New `bin/guess` md5
`8c1075a5928751333bbb111217a17c57` (gnu build of the same source as the coupled tree's
intel `4b46d356…`). The old binary is kept as `bin/guess.lpj_guess_awiesm3_20260705`. The
two runscripts still pin the coupled tree's binary; switch them to the setup's own `guess`
only after a smoke test of the gnu build.

**LUH3.** `lpj_guess.cmip.yaml` defaults to the 3-1-1 files in `land_use/regridded/`; the
coupled 11-series pins 3-1-2 in `land_use/TCO95/`. Both new spin-up runscripts pin 3-1-2 so
the state is made with what the coupled run reads.

**Leg length must be self-consistent.** esm_tools sizes the coupled leg from `nyear`, but
FESOM takes `run_length` from `restart_rate`. Keep `nyear`, `restart_rate` and
`restart_first` equal.

**A crashed leg can report success.** LPJ-GUESS ranks exiting 99 still left the job
`COMPLETED 0:0`, so esm_tools advanced the date file and resubmitted. `11X` produced three
legs in six minutes and zero years before it was caught. Check years-of-output, not exit
codes.

**Copies can truncate silently.** Staging `dist_1792` into `core3/` stopped at 159 of 3585
files with no error. Verify file counts after any large copy. (`/work` was 99 % then; 83 %
this evening.)

## What `11X` is for, once unblocked

Isolating whether corrected cavity depths fix the Weddell polynya diagnosed on 2026-08-26:
Oct-Dec `MLD2` at 90-60S is 758.6 m against WOA18's 148.8 m, a quarter of the band convects
past 1000 m and the deepest node reaches 4160 m, on a surface 0.20 psu too salty. That
convection mines heat from 700-2000 m, where the model is 0.358 K too cold against PHC3,
melts ice, drops albedo and raises absorbed shortwave. It is the mechanism behind the
positive TOA and its drift. It will settle if MLD2 falls toward 148.8, the >1000 m area
fraction falls from 25 %, and the upper-50 m salinity bias falls from +0.181 psu.
