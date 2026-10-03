# Handoff: winter sea ice too far south, Labrador Sea convection shutdown, and the 1200 s test

Written 2026-09-30 on albedo by the session that ran the albedo performance campaign
(`investigation_awiesm3_albedo_performance`). Everything below was measured; each claim has
its run, years and numbers. Run paths moved today: all project data now lives under
`/albedo/work/projects/p_awiesm3_cmip7/jstreffi/` (runtime, input, reval...), nothing is left
at the old top level.

## What to do

1. **Do not continue `PICAL_crunveg_ihf0` at `albsn` 0.82.** Its 2140s are entering the Arctic
   runaway (section 2). Jan decides; the default should be back to 0.80.
2. **Run `PICAL_crunveg_ts1200`** (section 5): the decade that decides whether the 1800 s FESOM
   time step or the new FESOM build shut down the Labrador Sea. Runscript is prepared and
   unsubmitted; one prerequisite (an OASIS restart set in `dist_1024` order) has to be made first.
3. **Do not resubmit the production runs** `PICAL_ccnice` (leg 2130) / `PICAL_crunveg` (leg 2120)
   before that test reports. Both runscripts are set up for the new FESOM build at 1800 s, i.e.
   the configuration that loses Labrador convection.

## 1. Runs involved (albedo)

Base `/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4/`.

| Experiment | Years done | FESOM build | FESOM step | Notes |
|---|---|---|---|---|
| `PICAL_crunveg` | 2100-2119 | leadclose 766020ec (libfesom 4c76f03d) | 1200 s | parent of everything below |
| `PICAL_crunveg_ihf0` | 2120-2149 | #1060+main `icehf_1e6cb53a` (de85c407); from 2140 `icehf_731d3030` (cbff2ecb, same numerics + shutdown fix) | 1800 s | branch of crunveg at 2119-12-31, OpenIFS cold start (discard 2120). `albsn` 0.80 for 2120s, **0.82 from 2130** |
| `PICAL_crunveg_ihf1` | 2120-2129 | icehf_1e6cb53a | 1800 s | as ihf0 + McPhee St\*u\* ocean-to-ice flux. Stopped: lowered Antarctic volume (Sep -1.4e3, Feb -0.4e3 km3) |
| `PICAL_crunveg_ihf0ctl` | 2130-2139 | icehf_1e6cb53a | 1800 s | branch of ihf0 at 2129-12-31, `albsn` 0.80: control for the 0.82 decade. OpenIFS cold start (discard 2130) |

Layout of all of them: 25 nodes, FESOM 768 x 1 at 96 ranks per node (`fesom.ranks_per_node: 96`,
esm_tools placement), OpenIFS 384 x 4, XIOS 16 x 4 two per OpenIFS node, rnfmap on an OpenIFS node.
~2.4 s per model day, a 10-year leg in 2 h 40 min.

## 2. Albedo: `albsn` 0.82 against the 0.80 control

Same decade, same parent state (ihf0 at 2129-12-31), 2131-2139, cold-start year discarded.
Volumes 1e3 km3, extents 1e6 km2 (FESOM `siextentn`/`sivoln`), `d/SE` = difference over the
standard error from interannual scatter (one realisation each, autocorrelation ignored).

| metric | ctl 0.80 | ihf0 0.82 | diff | d/SE | ihf0 0.82, 2141-49 |
|---|---|---|---|---|---|
| NH vol Apr | 32.46 | 35.15 | +2.68 | 5.4 | 35.38 |
| NH vol Sep | 11.77 | 13.78 | +2.02 | 4.5 | **16.07** |
| NH ext Mar | 18.11 | 18.59 | +0.48 | 2.0 | 17.40 |
| NH ext Sep | 6.90 | 6.95 | +0.05 | 0.3 | 7.66 |
| SH vol Feb | 1.31 | 1.44 | +0.13 | 0.4 | 1.55 |
| SH vol Sep | 13.47 | 13.37 | -0.10 | -0.2 | 13.93 |
| T2m 60-90N (C) | -9.04 | -9.67 | -0.64 | -4.1 | -9.74 |
| T2m global (C) | 13.87 | 13.73 | -0.13 | -3.9 | 13.73 |
| net TOA (W/m2) | 0.50 | 0.30 | | | 0.44 |

The Arctic takes the albedo step, the Antarctic does not move at all. Expected: the 0.75 -> 0.80
step gave the Antarctic about +1e3 km3 for +10 W/m2, so +0.02 is below one decade's detection
threshold (1.4e3 km3 SH Sep, 2-sigma from ihf0's scatter). **albsn cannot close the Antarctic deficit.**

Stop criteria written into ihf0's runscript (effective NH forcing = (albsn - 0.75) x
d(alb)/d(albsn)_NH x NH melt-season SWdn, from `scripts/analysis/ice_albedo_requirement.py` in
3-year windows; runaway reference: 0.83 grew 2.1 -> 7.3 W/m2 and momixoff ran away at 1.2 %
snow-free):

| window | run | albsn | NH snow-free | d/dalbsn NH | effective NH forcing |
|---|---|---|---|---|---|
| 2121-23 / 24-26 / 27-29 | ihf0 | 0.80 | 4.8 / 4.4 / 5.9 % | 0.238 / 0.249 / 0.260 | 2.05 / 2.19 / 2.31 |
| 2131-33 / 34-36 / 37-39 | ctl | 0.80 | 3.9 / 4.6 / 6.3 % | 0.257 / 0.246 / 0.244 | 2.2 / 2.1 / 2.1 |
| 2131-33 / 34-36 / 37-39 | ihf0 | 0.82 | 3.7 / **2.3** / 5.2 % | 0.342 / 0.403 / 0.273 | 4.24 / **5.07** / 3.40 |
| 2141-43 / 44-46 / 47-49 | ihf0 | 0.82 | **2.6 / 1.8 / 1.8 %** | 0.302 / **0.572** / 0.406 | 3.75 / **7.52** / 5.12 |

The 2140s cross both criteria in every window and reach 7.5 W/m2 at 1.8 % snow-free: the
amplifier that took 0.83 over the edge. 0.80 stays stable. PI caveat (memory
`pi-vs-present-day-ice-reference`): GIOMAS/PIOMAS are present-day, a PI Arctic must exceed them,
so some NH gain was wanted, but not a runaway.

## 3. Winter sea ice too far south (Jan noticed it on the ihf0 2135-2139 March map)

March maps for every run below are in
`/albedo/work/user/jstreffi/postprocessing/investigation_awiesm3_albedo_performance/figures/`
(tripyview, `hslice_np` a_ice settings, 5-year March means).

### 3a. Where the extra ice is (0.5-degree `a_ice.fesom.gr`, extent = cells >= 0.15, March, 1e6 km2)

| region | crunveg 2111-19 | ihf0 2121-29 | ihf0 2131-39 |
|---|---|---|---|
| NH total | 18.61 | 20.00 | 20.67 |
| Bering (50-66N, 163E-157W) | 1.26 | 1.82 | 2.16 |
| Labrador + Newfoundland (40-66N, 70-44W) | 1.30 | 1.72 | 2.07 |
| Okhotsk | 1.93 | 2.22 | 2.31 |
| GIN / Barents+Kara / Hudson / Baltic | 0.88 / 1.97 / 1.20 / 0.05 | 0.93 / 2.05 / 1.20 / 0.04 | 0.80 / 2.00 / 1.20 / 0.07 |
| south of 50N | 0.81 | 1.16 | 1.29 |
| southernmost ice, Bering box | 54.75N | 50.25N (box floor) | 50.25N |
| southernmost ice, Labrador box | 46.25N | 45.25N | 45.25N |

The central Arctic does not change; the growth is all in the marginal seas. FESOM's own March
NH extent: ccnice lineage 16.2-17.5 for 180 years, crunveg 16.5-16.6, **ihf0 2120s 17.88**,
2130s 18.58 (0.80 control 18.11). Note: the `.gr` files store lon as 0..360.

### 3b. History (levante + albedo), mesh-independent measure

Southernmost native-mesh node with March a_ice >= 0.15 in the Bering/N-Pacific box
(163E-157W, 40-66N); the Aleutians are at ~52N.

| run | mesh | years | median southmost | % years south of 52N |
|---|---|---|---|---|
| awicm3-v3.3.0 SPIN | CORE2 | 1350-1845 | 47.0N | 75 |
| awicm3-v3.3.0 PI | CORE2 | 1850-2014 | 47.3N | 65 |
| awicm3-develop BASE/ALB/L137 | CORE2 | 1990-2039 | 55-57N | 4-10 |
| awicm3-develop TUNE42PI_FES27 | CORE2 | 1350-1638 | 52.8N | 46 |
| awiesm3-develop momix_on/off_ssp_off | CORE2 | 1850-1859 | 45.7-47.3N | 80-90 |
| FESOM 2.6 standalone core3_spinup (JRA-type forcing) | old CORE3 (204875 nodes) | 1958-1997 | **56.4N** | **0** |
| awiesm3-v3.4 PI200 | CORE3 | 1390-1615 | 55.3N | 17 |
| PICAL_ccnice (levante) | CORE3 | 1940-2095 | 56.0N | 6 |
| PICAL_ccnice / crunveg (albedo) | CORE3 | 2100-2129 | 54-55N | 10-15 |
| **PICAL_crunveg_ihf0** | CORE3 | **2120-2129** | **50.1N** | **80** |
| ihf0, albsn 0.82 | CORE3 | 2131-2139 | 48.7N | 89 |

Pacific ice south of the Aleutians is an **old coupled-model bias**: strong already in the
first AWI-CM3 v3.3 spin-up (within its first 5 years), absent in the observation-forced ocean,
on both meshes. It was mostly gone in the PI200 / ccnice line and came back at 2120 with ihf0.
Other levante CORE3 runs (11X, 13A, 15F, 16C-E, the first decades of PICAL/momixoff) also had
it (decade-mean southmost 50-52N in the regridded scan). Cause on the Pacific side is not
investigated: no deep convection there, so section 4 does not apply; candidates are winter
winds / cold air over the Bering Sea or ice export. Labrador reaches ~45N in nearly every run
including the forced ocean (45.8N), so the southern Labrador ice *edge* is closer to a baseline
feature; what changed there is the interior (section 4).

## 4. The Labrador Sea stopped convecting

The TUNE42PI_FES27 map (AWI-CM3, CORE2, KPP, K_GM 1000) has an open Labrador Sea interior in March;
ihf0 has it ice-covered. March means, last 5 years, native mesh (`MLD2` absolute):

| region | run | ice | MLD mean | % nodes > 1000 m | SSS |
|---|---|---|---|---|---|
| Labrador (56-62N, 60-50W) | TUNE42PI_FES27 1635-39 | 0.40 | 1613 m | 54 | 34.51 |
| | ccnice levante 2095-99 | 0.39 | 835 m | 36 | 34.17 |
| | crunveg 2115-19 | 0.38 | 1076 m | 44 | 34.29 |
| | **ihf0 2135-39** | **0.62** | **100 m** | **0** | **33.88** |
| Irminger (58-63N, 40-30W) | TUNE42 / ccnice / crunveg / ihf0 | 0.02-0.03 | 2062 / 655 / 293 / 220 m | 78 / 26 / 0 / 0 | 35.07 / 35.00 / 34.81 / 34.41 |
| Greenland Sea (72-77N, 10W-5E) | TUNE42 / ccnice / crunveg / ihf0 | 0.67 / 0.30 / 0.16 / 0.18 | 663 / 1838 / 2344 / 1217 m | 27 / 89 / 90 / 67 | 34.48 / 35.01 / 35.22 / 35.12 |

Labrador interior, March, year by year:

| run | years | MLD | % > 1000 m | SSS | ice |
|---|---|---|---|---|---|
| crunveg (leadclose, 1200 s) | 2100-2119 | 290-1915 m, typically ~1 km | 0-63, mostly 30-60 | 33.95-34.74 | 0.09-0.56 |
| ihf0 (new build, 1800 s) | 2120-2127 | 254-639 m | 5-22 | 33.91-34.51 | 0.26-0.53 |
| ihf0 | **2128-2139** | **56-160 m** | **0 every year** | 33.69-34.00 | 0.20-0.88 |
| ihf0ctl (albsn 0.80) | 2130-2139 | 65-322 m | 0-8 | 33.53-34.26 | 0.21-0.92 |
| ihf0 (albsn 0.82) | 2140-2149 | 52-124 m | 0 | 33.44-33.97 | 0.27-0.87 |

- Weakened **immediately at the 2120 branch**, collapsed in **2128**, has not recovered in 22 years.
- **Not the albedo**: the collapse precedes the 0.82 step (2130), and the 0.80 control is shut down too.
- It locks in: capped surface ~0.5 fresher than crunveg, no convective supply of heat/salt,
  basin freezes; melt of the extra ice freshens it further. The Irminger Sea had already stopped
  in crunveg and keeps freshening; the Greenland Sea still convects.
- At the 2120 branch **exactly two things changed** (namelist diff of crunveg's 2110 leg against
  ihf0's 2120 leg): FESOM `step_per_day` 72 -> 48 (1200 -> 1800 s) and the FESOM build
  (leadclose 766020ec -> #1060+main 1e6cb53a with use_ustar_ocnice off). Everything else in
  every FESOM namelist is identical; lead closing (`h0cur = clamp(hold, h0min, h0max)`,
  0.5/1.5) is identical code in both builds.
- The 2-month tests from the same 2100-01-01 restart cannot separate them: all convect to
  1.6-1.9 km in their first winter (perfE/G 1200 s old build 1704 m; perfF/H 1800 s old build
  1772 m; perfI 1800 s new build 1636 m; nolockB1/B4, identical setup, 1887 / 1735 m). The
  shutdown builds over years.
- Differences TUNE42 (convects) against the AWI-ESM3 CORE3 line, for later: CORE2 vs CORE3+cavities,
  zstar vs linfs, KPP vs cvmix TKE+IDEMIX, K_GM_max 1000 vs 2500 / Redi_Kmax 0 vs 1000, AWI-CM3
  vs AWI-ESM3. Ready-made one-at-a-time sensitivity runs on levante:
  `PICAL_cvKPP / cvTKE / cvTKEIDEMIX` (2020-2039), `PI200_gmR1500 / gmR2500` (1580-1599),
  `16C` (K_GM 1500) / `16E` (2500) (1350-1389). Not yet analysed for Labrador convection.

## 5. The test: `PICAL_crunveg_ts1200`

**Question**: does the 1800 s FESOM step, or the new FESOM build, shut down the Labrador Sea?

**Design**: identical to ihf0's first decade (branch of `PICAL_crunveg` at 2119-12-31, new
FESOM build, constant gamma_t, albsn 0.80, OpenIFS cold start) except `fesom.time_step: 1200`.
One decade, 2120-2129; ihf0 collapsed within 8 years. Compare with ihf0 2120s (same build,
1800 s) and crunveg 2100s/2110s (old build, 1200 s).

**Decision rule** (Labrador interior, March, 2121-2129, `scripts/analysis/seaice_20260930/labrador_yearly.py`):
- keeps convecting like crunveg (MLD ~1 km, 30-60 % of nodes > 1000 m most years): **time step**
  is the cause -> production back to 1200 s with the layout below.
- weakens/collapses like ihf0 (MLD a few hundred m, then < 150 m): **the FESOM build** is the cause
  -> bisect the build (main vs leadclose) next; the time step is free to stay at 1800 s.
- ambiguous (partial weakening): run the complementary arm (old build at 1800 s, below).

**Runscript** (prepared, not submitted):
`/albedo/home/jstreffi/esm_tools/runscripts/awiesm3/develop/awiesm3-develop-albedo-TCO95L91-CORE3_PICAL_crunveg_ts1200.yaml`
It is `PICAL_crunveg_ihf1.yaml` (first decade) with six changes: time_step 1200, FESOM nproc 1024,
ranks_per_node 64, `use_ustar_ocnice: .false.`, FESOM build `fesom-bin/icehf_731d3030`
(libfesom cbff2ecb: ihf0's code + the profiler-after-MPI_Finalize fix, verified: ihf0's 2140
leg is the first new-build leg whose srun step ends COMPLETED), and the OASIS
`ini_restart_dir` below.

**Layout for 1200 s (33 nodes)**, estimated from measured LUCIA FESOM compute (2-month tests):
new build at 1800 s: 768 @ 64/node 135-138 s, 768 @ 96/node 155-158 s; 1800 -> 1200 s factor
1.37-1.39 (perfE/G vs perfF/H, OpenIFS unchanged); 768 -> 1024 ranks at equal density 0.895
(perfJ vs perfM); OpenIFS compute 166-177 s.

| FESOM at 1200 s | nodes | est. FESOM compute | vs OpenIFS ~170 s |
|---|---|---|---|
| 768 @ 96/node | 25 | ~215 s | FESOM limits +26 % |
| 768 @ 64/node | 29 | ~189 s | FESOM limits +11 % |
| 1024 @ 96/node | 28 | ~193 s | FESOM limits +13 % |
| **1024 @ 64/node (chosen)** | **33** | **~169 s** | **balanced** |

Expect roughly 2.5 s per model day (a 10-year leg ~2 h 45 min). Unmeasured: 1024 single-threaded
spread and the 1200 s factor on the new build; if a speed check is wanted first, a 2-month TEST
with `oasis3mct.use_lucia: true` and `bench/lucia_steps.py` from the performance repo gives it.
The science run itself keeps LUCIA off.

**Prerequisite: an OASIS restart set in `dist_1024` order from crunveg's 2119-12-31 state.**
The runscript points at `/albedo/work/projects/p_awiesm3_cmip7/jstreffi/input/oasis3mct_reorder/PICAL_crunveg_dist1024_21191231/`,
which does not exist yet. Source: `.../oasis3mct_reorder/in_PICAL_crunveg_dist1792/` (crunveg's
original end-2119 set, levante dist_1792 order). Recipe (memory `oasis-reorder-tool-albedo`,
read the tool's README in full first):
- tool: `/albedo/work/user/jstreffi/software/oasis_reorder_tool` (built); run under
  `srun -p smp --qos=30min -A paleodyn.paleodyn -n 11 -c 4` with `ulimit -s unlimited`;
  source partition = levante dist_1792 (the project mesh dir `input/fesom2/core3/dist_1792` IS
  levante's), target = `input/fesom2/core3/dist_1024`. The earlier
  `TCO95-CORE3_dist1024_from_levante1792` set was made exactly this way (logs `reorder*.log` there).
- the tool skips grids/masks/areas, but the atm->feom couplings are `CONSERV` and read FESOM
  areas/masks: copy them in as `rmp_zz_{grids,masks,areas}.nc`, reorder, rename back.
- check the order: grids.nc feom lon/lat against `nod2d.out` in the concatenated `my_list` order
  of dist_1024 (1609 polar nodes differ < 0.006 deg by rotation roundoff; a wrong order is off
  by up to 180 deg everywhere).
- add dated aliases so oasis.py's glob `f"{ini_restart_dir}{restart_file}*{YYYYMMDD}"` finds them:
  `ln -s rstas.nc rstas.nc_21100101-21191231` and the same for `rstos.nc`, `vegin.nc`.
- FESOM restarts come from crunveg's own `restart/fesom/` (partition-independent netCDF), LPJ-GUESS
  from `lpjg_state_2120`; both already set as in ihf1.

**Submit**: `esm_runscripts <runscript> -e PICAL_crunveg_ts1200 --no-motd` (new experiment, no
prompt; `use_venv: True` builds `<exp>/.venv_esmtools` from GitHub `feat/awiesm3-v3.4-co2`,
322f0ed1f, ~3.5 min; it has the `ranks_per_node` placement). After prepare, check in
`run_21200101-21291231/work`: `namelist.config` step_per_day 72, `rank_placement` 1489 lines
with at most 64 FESOM ranks per node, `hostfile_srun`, `#SBATCH --nodes=33`, libfesom md5
cbff2ecb, `rstas.nc` md5 equal to the reordered set, `fesom.2120-01-01.*.restart` staged.
Discard 2120 in the evaluation.

**Complementary arm, only if needed**: old build at 1800 s = crunveg branch with FESOM from
`PICAL_crunveg/bin/fesom` (libfesom 4c76f03d, leadclose 766020ec) and `time_step: 1800`,
768 @ 96/node, OASIS from crunveg's `restart/oasis3mct/` (dist_768 order, the ihf runs used
it). With both arms the 2x2 (build x step) is complete.

## 6. Things that cost time today (do not repeat)

- **Project path moved**: everything is under `/albedo/work/projects/p_awiesm3_cmip7/jstreffi/`;
  runscripts, finished configs, symlinks and the esm_tools venvs were rewritten/rebuilt.
- **Experiments stage FESOM from their own `<exp>/bin/` copy** taken at creation: changing
  `bin_sources` later does nothing without `--update-filetypes bin` (verify md5 in `work/lib/fesom/`).
- **venv is frozen**: esm_tools changes need a push to `feat/awiesm3-v3.4-co2` and the venv moved aside.
- **Resubmitting an existing experiment with a changed runscript**: use `-U` (else it prompts).
- **Branch legs cold-start OpenIFS**: discard the first year; continuation legs are warm.
- **`vegin.nc` in a branch's first leg** comes from the pool (LPJ-GUESS input overrides the restart copy); affects only the first daily exchange.
- **FESOM builds without 731d3030** end every leg with `MPI_Comm_f2c() called after MPI_FINALIZE`
  (profiler report after oasis_terminate; FESOM_PROFILING defaults ON in main). Harmless for
  output/restarts, step shows CANCELLED. Fix branch pushed: FESOM `fix/profiler-report-before-mpi-finalize` (no PR yet).
- **XIOS**: largest server 55.9 GB per 10-year leg; keep 2 per OpenIFS node.
- **Reading albedo XIOS output** needs `HDF5_PLUGIN_PATH=/albedo/work/user/jstreffi/model_codes/awiesm3-develop/xios/hdf5-blosc/plugin` (Blosc filter 32001); threaded netCDF opens segfault without it.
- **tripyview**: env `esm-tools_auto_tripyview` (albedo `~/.conda/envs`, levante too) repaired to
  tripyview 1aeabd7 (`/albedo/work/user/jstreffi/software/tripyview`; on levante `~/tripyview_1aeabd7`).
  Needs whole-file chunks, `data[["a_ice"]]` (XIOS bounds variables), and years < 1678 break it
  (see `tripy_pi200.py`). The esm_tools Tripy_Fesom subjob has never run on albedo:
  fesom-2.7/2.8.yaml set plain `scripts_files/sources/targets` for the XIOS rename scripts and
  override fesom.tripyview.yaml's; they should be `add_scripts_*`.

## 7. Scripts and data

`scripts/analysis/seaice_20260930/` in this repo:
- `labrador_yearly.py`, `convection_boxes.py`, `labrador_tests.py`: convection diagnostics (edit run/years inside or pass args).
- `ice_regions.py`: regional 0.5-degree extents; `ice_south_native.py` / `ice_south_levante.py`: ice-edge scans.
- `tripy_aice_run.py`, `tripy_aice_np.py`, `tripy_pi200.py`: the March maps.
- `ice_south_levante_all.txt`, `ice_native_all.txt`, `native_albedo.txt`: raw scan output.
Evaluation outputs of the albedo decade (csv, `alb_req_*.txt`) are in
`/albedo/work/user/jstreffi/postprocessing/investigation_awiesm3_albedo_performance/.eval_tmp/`.

## 8. Addendum 2026-09-30 (evening): the test is now a full 2x2 in the same decade

**Why the decision rule in section 5 was not enough.** Labrador convection is intermittent in
the old configuration too. `PICAL_ccnice` (old build, 1200 s, albedo) has single-winter
collapses to ~100 m in 2108, 2112, 2117, 2121, 2126 and 2129, more frequent in the 2120s,
each recovering within 1-2 winters; crunveg was sturdier (weak only 2105, 2117). So crunveg's
2100s are not a control for 2120-2129 in a drifting system, and a single weak winter is not
the signal. What separates the new-build 1800 s runs is that they **never recover**: ihf0
from 2128 (to 2149), ihf1 (the McPhee arm, second realisation) from 2127.

**The extra Labrador ice is a consequence, not the trigger.** FESOM's own tendencies, Nov-Apr,
Labrador interior: crunveg imports 86 cm and melts 72 cm; ihf0 imports 70 and melts 53, ihf1
82 and 54. Less import, less melt: the ice stays because convection no longer brings heat up.
The Bering extra ice is different, locally frozen (thdgrice 28 -> 35 cm).

**The build diff is small.** 766020ec -> 53e095ec touches 19 files; at one FESOM thread the
ice rheology, advection and IDEMIX changes are refactors (lock removal, node-owned gathers,
EVP halo overlap 434e2fcd, one OpenMP region). `8d2cd552 ifs: restore legacy sea-ice
coupling` changes only the `__ifsinterface` path, not the `__oifs` one we build.

**Runs (all branches of PICAL_crunveg at 2119-12-31, OpenIFS cold start, discard 2120):**

| | 1200 s | 1800 s |
|---|---|---|
| new build (icehf 731d3030 / 1e6cb53a, same numerics) | `PICAL_crunveg_ts1200` (1024 @ 64, 33 nodes) | `PICAL_crunveg_ihf0` (768 @ 96, 25 nodes), done |
| old build (leadclose 766020ec, libfesom 4c76f03d) | `PICAL_crunveg_ob1200` (as ts1200) | `PICAL_crunveg_ob1800` (as ihf0) |

Each old-build arm differs from its new-build partner only in the FESOM binary
(`fesom-bin/leadclose_766020ec/`, copied from crunveg's own per-leg library) and in having no
`use_ustar_ocnice` line: the old build's `/ice_therm/` does not declare it and the read
swallows errors, so the line would silently reset the whole group. Staged namelists checked
against the partner: that line is the only difference.

**OASIS for 1024**: `input/oasis3mct_reorder/PICAL_crunveg_dist1024_21191231/`, reordered from
`in_PICAL_crunveg_dist1792` with `reorder_oasis` (11 files, grids/masks/areas via `rmp_zz_*`),
byte-identical weights/grids/masks/areas to the earlier dist1024 set, grids.nc feom coordinates
match nod2d.out in dist_1024 my_list order on all 220509 nodes (max 0.006 deg); negative
control (1792 source against dist_1024) fails on every node. Check script:
`scripts/model/check_oasis_order.py`.

**Score** with `scripts/analysis/seaice_20260930/labrador_2x2.py`: collapsed winters, the longest
consecutive run, and whether the arm ends convecting. If the 1200/1800 or old/new contrast is
still ambiguous after 2121-2129, extend the ambiguous cells to 2139 (ihf0 and ihf1 locked in
only in year 8-9).

**Underneath all of it**: the subpolar gyre has been freshening across the whole lineage
(Irminger SSS 35.07 TUNE42, 35.00 ccnice, 34.81 crunveg, 34.41 ihf0; Irminger convection
already off in crunveg), consistent with levers adopted in rounds 36-39 (GM 2500 weakens the
AMOC, TKE+IDEMIX grows Arctic ice). Whichever factor wins the 2x2, the Labrador Sea stays one
perturbation from switching off until that drift is dealt with.
