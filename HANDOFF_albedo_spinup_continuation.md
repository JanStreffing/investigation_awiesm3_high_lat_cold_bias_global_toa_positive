# Handoff: run the AWI-ESM3 pre-industrial spin-up on albedo

Rewritten 2026-09-25 18:00 on levante, for whoever picks this up on albedo.

`PICAL_ccnice` is the AWI-ESM3 v3.5 pre-industrial spin-up line: OpenIFS 48r1 TCO95L91 +
FESOM2 CORE3 with Antarctic ice-shelf cavities + LPJ-GUESS + OASIS + XIOS. It has run its
own years 1940 to 2100 on levante. The reason for moving to albedo is queue time, not
anything wrong with the run: a 47-node leg on levante has waited anywhere from 30 seconds
to 10 hours to start, and 50 model years is five such legs.

Read **Three things that will kill your first run** before submitting anything. All three
are one-line problems and all three are fatal.

## Where the run actually is, and the branch point question

Complete through 2099. `fesom.2100-01-01` restarts written. As of this writing two 50-year
continuations are **queued but not started** on levante, both 2100 to 2150, both from the
same 2099-12-31 restarts:

| run | land it starts from | purpose |
|---|---|---|
| `PICAL_ccnice` | its own, carried forward | the continuation |
| `PICAL_crunveg` | `lpjg_state_3850`, the offline CRUNCEP spin-up end state | gives the canopy work a forest to act on |

Both were submitted on 2026-09-25 with estimated starts around 02:30 the next morning, and
**both are expected to be cancelled** in favour of moving the work here. Do not assume they
produced anything. Check before relying on a restart later than 2100:

```
ls /work/bb1469/a270092/runtime/awiesm3-v3.4/PICAL_ccnice/restart/fesom/
ls /work/bb1469/a270092/runtime/awiesm3-v3.4/PICAL_crunveg/restart/fesom/ 2>/dev/null
```

**So the branch point is 2100 unless those directories say otherwise.** That is the only
point where a complete, consistent restart set exists, and if the levante runs are cancelled
there is no duplication to weigh: albedo simply becomes the line.

What does not go away if they are cancelled is the question they were asked to answer, since
both are cheap to reproduce here and both are worth having. `PICAL_crunveg` in particular is
the only test of whether the CRUNCEP state loads on `TCO95-land` at all, which is a
prerequisite for any land experiment on this configuration. Its runscript is
`awiesm3-develop-levante-TCO95L91-CORE3_PICAL_crunveg.yaml` and differs from the
continuation in exactly one thing, the LPJ-GUESS state it starts from.

## The state you are inheriting

Measured over 2060-2099, the whole `cvmix` era. This is **not** a converged control and any
evaluation should say so.

| metric | 2060-64 | 2095-99 | trend/decade | 95% CI |
|---|---|---|---|---|
| T2m global | 14.439 | 13.645 | **-0.214** | 0.025 |
| T2m 60-90N | -8.143 | -9.842 | **-0.511** | 0.118 |
| NH ice extent, March | 16.262 | 18.188 | +0.522 | 0.103 |
| net TOA | -0.775 | +0.236 | +0.272 | 0.079 |

Every one is outside its confidence interval. The run cools 0.21 K/decade globally while the
Arctic grows ice. Net TOA over the last two decades is +0.167 W/m2 and total ocean heat
content is flat (+0.04 ZJ/decade), so the cooling is not an energy leak: the upper 700 m
loses heat to below 700 m at about 0.32 W/m2 equivalent, almost exactly balanced.

Against GIOMAS 1989-1999 (use `OBS_Y0=1989 OBS_Y1=1999`; the early record is the fairer
anchor for a PI run):

| 10^3 km3 | model 2090-99 | GIOMAS 1989-99 | bias |
|---|---|---|---|
| NH April volume | 34.56 | 31.19 | +3.37 |
| SH September volume | 12.33 | 19.89 | **-7.56** |
| SH February volume | 0.73 | 2.07 | **-1.34** |

The Antarctic sits at about 62 % of observed winter volume and has been the campaign's open
problem for months. Drake Passage transport is 98.9 Sv against Cunningham 2003's 134 +- 11.2,
and it has been drifting down about 7 Sv over 2040-2099. GM tuning cost 16.7 Sv of that:
`PI200` 110.0, `gmR1500` 103.2, `gmR2500` 93.3 Sv on identical years, monotonic in
`K_GM_max`. That was not scored when GM 2500 was adopted for the CDW improvement.

CMPI went 0.844 (1980-89) to 0.903 (2090-99), and the degradation is **tropical**: nino34
+0.155, tropics +0.134, against arctic +0.003 and antarctic +0.008.

**The boreal forest is gone.** In the Siberian box 55-75N/60-180E at 2099: BNS 0.001,
TREEFPC 0.050, GRASSFPC 0.216, grass beats tree on 84.9 % of cells, against 0.227 / 0.415 /
29.2 % in the CRUNCEP spin-up state it descends from. It is lost to competition with C3
grass, not to a climate gate: every BNS establishment gate stays open while it disappears,
the only temperature gate is a warm one 34 K away, and GDD5 clears on 91 % of cells. Do not
answer a question about the boreal forest through thresholds.

## The software stack, which you should not improvise

`esm_tools`: branch `feat/awiesm3-v3.4-co2`, at **`588cf7af7`**. It carries the v3.5 CORE3
defaults this run depends on: `cavity_gamma_scale: 0.6`, sea-ice albedos (`albsn` 0.80,
`albsnm` 0.65, `albi` 0.70, `albim` 0.68, `albpnd` 0.20), `use_momix: false`,
`RCCNSEAICE: 15.0`, `RCL_INPPMIN: 70000.0`, `mix_scheme: cvmix_TKE+cvmix_IDEMIX` and the
IDEMIX forcing paths. The runscript sets `install_esm_tools_branch` to this branch.

| component | branch | note |
|---|---|---|
| oifs-48r1 | `awiesm3-develop` | fine as installed |
| **lpj_guess** | **`lpj_guess_awiesm3`** @ `1325440` | changed 2026-09-25, see below |
| xios | `main` | fine |
| rnfmap | `enthalpy_on_ocean_side` | fine |
| oasis3mct | `local_combined_fixes` (SMHI) | fine |
| **fesom** | see below | **does not work as installed** |

**FESOM.** `esm_master install-awiesm3-develop` resolves `fesom-2.7-main` to
`feat/cavity-gamma-scale`, which is 53 commits behind `main` and **will not run this
experiment**: it has `cavity_gamma_scale` and `use_atm_ice_tskin` but neither IDEMIX fix
(three `MPI_DOUBLE_PRECISION` reductions remain, no `nf90_fill_float` fallback).

You need `main` + PR #1060 + PR #1087. Verified to merge with one trivial conflict in
`src/gen_modules_config.F90`, where `main` and #1060 each appended a name to the same
`namelist /run_config/` continuation line; keep both. **Check first whether they have
merged** — as of 2026-09-25 both are open with Suvarchal Cheedela and Patrick Scholz as
reviewers. If they have, point `fesom-2.7.yaml`'s `2.7-main` at `main` and stop there;
`feat/cavity-gamma-scale` has outlived its purpose either way.

**LPJ-GUESS changed on 2026-09-25** and this is the newest thing here. `lpj_guess_awiesm3`
now carries the AWI land-grid line merged with Laszlo Hajdu's canopy work: stem area index
(via `ltos`), shrubs routed to low IFS vegetation, `iftreefracca 1`, his nitrogen roundoff
and land-use transfer fix, alongside the Raupach roughness sent to OpenIFS as `GUE_Z0HV`,
the LUH3 crop fractions and the coastline land mask. Built as
`guess.4112_awiesm3canopy_1325440_20260925`.

The grid is now chosen by configuration rather than hard-coded: the code reads
`grid_product` and `grid_name` from `lpjg_steps.yaml`, which `lpjg_steps.yaml.j2` already
writes as `coupled_land` and `TCO95-land`. That needs no work from you; it is already
consistent.

## Three things that will kill your first run

**1. `h0min` in the runscript.** It sets `h0min: 0.5`, which is the code default, so the line
does nothing. But `h0min` is only in the `namelist /ice_therm/` group on the abandoned
`feat/leadclose-h0min` branch. On `main` and on the integration branch it is a plain local in
`ice_thermo_cpl.F90`, and FESOM aborts at the namelist read. **Delete the line.**

**2. The IDEMIX forcing files.** `awiesm3.yaml` reads both from the FESOM pool:

```
${fesom.pool_dir}/forcing/idemix/fourier_smooth_2005_cfsr_inert_rgrid.nc
${fesom.pool_dir}/forcing/idemix/forcing_idemix_final_bin/FIN_STORMTIDE2_M2_plus_NYCANDER_CnoM2_bin_0.40deg.nc
```

with `fesom.pool_dir` = `${general.pool_dir}/fesom2/`, and `general.pool_dir` on albedo is
`/albedo/pool/`. Both files must exist in that layout before the first submit. The paths the
stock FESOM 2.7 namelist ships are useless: its surface file is a DKRZ `/pool/data/` path and
its bottom file an `/albedo/pool/FESOM/fesom2.0/` path, so neither machine resolves both.

**3. `ltos` must be in the instruction file the run actually reads.** That file is
`namelists/lpj_guess/global.ins`, reached via `run_coupled_4_1_2.ins`. It is **not**
`ecearth.ins.j2`, which looks like the right file and is not used by this setup. `global.ins`
needs `ltos 0.1` in `group "common"` and `ltos 0.05` in `group "grass"`; without them the
canopy code is compiled and never exercised, silently. It is correct in `588cf7af7`; verify
it survived into the staged copy in your run directory.

## The gate: check this before you let a leg run

Within the first minute of the compute log:

```
--> IDEMIX total srf. energy Etot_srf =  0.3211645      TW
--> IDEMIX total tidal energy Etot_bot =  0.9135067     TW
--> IDEMIX Etot_bot after normalizing =  0.9179869      TW
```

0.9179869 is the forcing file's own `Etot_in_TW_bin` attribute, so this is a real check, not
a self-consistency one. `Etot_srf = 6.05e-04` with `Etot_bot = -1.54e16` means a FESOM
without the fixes. **Kill the job.** That signature NaN'd the ocean in 100 seconds on
2026-09-24 and cost a 47-node slot.

## Facts you should not have to rediscover

- **The staged library can be a symlink, and `md5sum` follows it.**
  `<exp>/bin/fesom/lib/fesom/libfesom.so` is what the run loads, not the `bin_sources` path.
  esm_tools keeps per-leg snapshots `libfesom.so_<leg>` and points `libfesom.so` at one. On
  2026-09-24 it pointed at a pre-fix snapshot while every checksum looked right. Use `ls -l`
  before `md5sum`. The `fesom` executable is a driver stub, byte-identical across all these
  builds, so a matching binary checksum proves nothing. Only the `.so` distinguishes them.
- **Chained legs do not read runscript edits.** `<exp>/scripts/<name>.yaml` is a symlink into
  a per-leg archive. Edit the source and submit with `-U`.
- **The date file advances past a crashed leg.** Count years in `outdata/` instead.
- **Preserve failed run directories.** Move them aside with a descriptive name, never `rm -rf`.
- **Never run two `esm_master` in one `model_codes` directory.**
- **OASIS `rmp_*.nc` weights belong to one FESOM partition.** Changing `fesom.nproc` is safe
  because the lookup is keyed on it; an `oasis3mct.ini_restart_dir` carrying `rmp` files will
  shadow that lookup silently.
- **Contained-run:** `general.use_venv: True` equals `--contained-run`, but `--open-run` in
  `original_command` overrides the yaml and is self-perpetuating. Break the chain and
  resubmit fresh.
- Do not submit `esm_runscripts` from inside a Slurm job; inherited `SLURM_MEM_PER_NODE`
  kills `srun`.
- **tripyview** is installed in `~/.conda/envs/esm-tools_auto_tripyview` but its `egg-link`
  points at a path that no longer exists. Run with
  `PYTHONPATH=/work/ab0246/a270092/software/tripyview`.

## Running it

`nyear: 10` and `restart_rate: 10` must match or the leg cannot finish. On levante the run
used 47 nodes at about 77 SYPD under TKE+IDEMIX (3:18:32 for 2080-89, 2:58:12 for 2090-99;
about 86 SYPD for `cvmix_TKE` alone and 75 for KPP). Albedo's node layout and core count
differ: treat those as a levante baseline and measure your own.

What has to travel, from `/work/bb1469/a270092/runtime/awiesm3-v3.4/PICAL_ccnice/`:

- `restart/fesom/fesom.<year>-01-01.{ice,oce}.restart`
- `restart/oifs/<year-1>1231/` — OIFS restarts are directories named by the **last day** of
  the leg, so 2100 means `20991231`
- `restart/lpj_guess/lpjg_state_<year>`
- `restart/oasis3mct/` for the same date
- the LPJ-GUESS **input** state under
  `/work/ab0246/a270092/input/lpj-guess/restart/TCO95-CORE3/`, which is separate from the
  run's own restart tree and easy to miss
- the runscript,
  `~/esm_tools/runscripts/awiesm3/develop/awiesm3-develop-levante-TCO95L91-CORE3_PICAL_ccnice.yaml`,
  which is deliberately **not** tracked in git

Evaluation lives in the campaign repo on levante,
`/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive`:
`scripts/analysis/ccnice_regional_drift.py`, `ice_volume_series.py`,
`giomas_volume_targets.py`, `drake_passage_transport.py`, `coupled_annual_store.py`, and
`scripts/figures/`. The full `reval` tool is
`/work/ab0246/a270092/software/release_evaluation_tool2`; copy
`configs/awiesm3_v35_pical_ccnice_2090s.py` and change `model_version`,
`pi_ctrl_start`/`end`, `spinup_end` and `out_path`.

## Things not to do

- **`surf_relax_s` stays 0.** Sea surface salinity restoring is excluded as a lever, even for
  Weddell convection. Settled repeatedly.
- **No hemispheric tuning parameters.** `albsn` and `albsnm` are global. One acts mainly in
  the south and the other mainly in the north because of measured surface state
  (`d(alb)/d(albsn)` is 0.908 in SH DJF against 0.149 in NH JJA), not a hemispheric switch.
- **Price albedo changes in W/m2 first.** The Arctic runs away above about +5.5 W/m2, the
  Antarctic is stable to at least +10.4. The sensitivities are state-dependent and have moved
  by a factor of three as the Arctic cooled: re-measure rather than reusing a number.
- **Do not open pull or merge requests unprompted** in any component repository. Branches,
  commits and pushes are fine. Jan controls what goes upstream.

## Open questions

**Does the canopy work do anything here?** It acts on the tree side of the tile and this run
has almost no trees. Laszlo measured +0.46 K global T2m from it, concentrated in the boreal
zone, in a run that has a forest. `PICAL_crunveg` is the test; until it reports, treat the
+0.46 K as not transferable to this configuration.

**Will the CRUNCEP state load on `TCO95-land`?** Unresolved at the time of writing. That state
was written on an 11538-cell gridlist from `L096.msk` and this configuration uses
`TCO95-land`, so LPJ-GUESS meets cells the spin-up never simulated. The old
`HANDOFF_lpjg_respin_core3.md` recorded that as fatal with peatland on; the code has changed
substantially since and the claim is untested. `PICAL_crunveg`'s first leg settles it.

**NH April volume** decides whether `cvmix_TKE+cvmix_IDEMIX` keeps its place. It was
+1.145 +- 0.968 per decade over 2080-2099 and sits 3.37 above the early GIOMAS record.

**The upper-ocean heat drain** is the mechanism behind the cold drift and nobody has tried to
stop it. Whether it is the mixing scheme, the Southern Ocean freshwater regime, or both is
open.

**The Southern Ocean deficit is not a vertical-mixing problem.** Two regimes are
characterised (salty/venting and fresh/warm CDW); `spp`, `h0`, MLE, stronger GM and
basal-melt 0.6 were all tested in 2026-09 without fixing it in 20 years.
