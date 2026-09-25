# Handoff: continue the PICAL_ccnice spin-up on albedo

Written 2026-09-25 on levante, for whoever picks this up on albedo.

`PICAL_ccnice` is the AWI-ESM3 v3.5 pre-industrial spin-up line: OpenIFS 48r1 TCO95L91 +
FESOM2 CORE3 with Antarctic ice-shelf cavities + LPJ-GUESS + OASIS + XIOS. It ran its own
years 1940 to 2100 on levante and stopped because it reached the `final_date` it was given,
`2100-01-01`. Nothing is broken. Your job is to stand the same experiment up on albedo and
carry it forward from 2100.

Read "Two things that will kill your first run" before you submit anything. Both are one-line
problems and both are fatal.

## The state you are inheriting

Measured on levante over 2060-2099, which is the whole `cvmix` era. This is **not** a
converged control, and any evaluation you write should say so.

| metric | 2060-64 | 2095-99 | trend/decade | 95% CI |
|---|---|---|---|---|
| T2m global | 14.439 | 13.645 | **-0.214** | 0.025 |
| T2m 60-90N | -8.143 | -9.842 | **-0.511** | 0.118 |
| NH ice extent, March | 16.262 | 18.188 | +0.522 | 0.103 |
| net TOA | -0.775 | +0.236 | +0.272 | 0.079 |

Every one of those is outside its confidence interval. The run is drifting cold at about
0.21 K/decade globally while the Arctic grows ice. Net TOA is slightly positive over the
last two decades (+0.167 W/m2) and the total ocean heat content is flat (+0.04 ZJ/decade),
so the cooling is not an energy leak. It is a redistribution: the upper 700 m loses heat to
below 700 m at about 0.32 W/m2 equivalent, almost exactly balanced.

Against GIOMAS 1989-1999 (the early record is the fairer anchor for a PI run; use
`OBS_Y0=1989 OBS_Y1=1999` with `scripts/analysis/giomas_volume_targets.py`):

| 10^3 km3 | model 2090-99 | GIOMAS 1989-99 | bias |
|---|---|---|---|
| NH April volume | 34.56 | 31.19 | +3.37 |
| NH September volume | 13.97 | 11.96 | +2.01 |
| SH September volume | 12.33 | 19.89 | **-7.56** |
| SH February volume | 0.73 | 2.07 | **-1.34** |

The Antarctic is at roughly 62% of observed winter volume and has been the campaign's open
problem for months. The Arctic is over-grown and still growing. CMPI over the same window
went 0.844 (1980-89) to 0.903 (2090-99), and the degradation is **tropical**, not polar:
nino34 +0.155 and tropics +0.134, against arctic +0.003 and antarctic +0.008.

Vertical mixing is `cvmix_TKE + cvmix_IDEMIX` (`mix_scheme` 56) from 2080. It was adopted
because it grew the Antarctic summer pack in a 20-year test arm, and on the production line
it did not reproduce that. It did repair an energy deficit that `cvmix_TKE` alone had opened
(net TOA -0.54 and -0.39 W/m2 in the 2060s and 2070s). Cavity basal melt under it is
1062 Gt/yr, against Rignot 2013 ~1500 and Adusumilli 2020 ~1100.

## What has to travel from levante

Experiment root on levante: `/work/bb1469/a270092/runtime/awiesm3-v3.4/PICAL_ccnice/`

- `restart/fesom/fesom.2100-01-01.ice.restart` and `.oce.restart`
- `restart/oifs/20991231/` (OIFS restarts are directories named by the last day of the leg)
- `restart/lpj_guess/lpjg_state_2100`
- `restart/oasis3mct/` for the same date
- The LPJ-GUESS **input** state the runscript points `ini_restart_dir` at, under
  `/work/ab0246/a270092/input/lpj-guess/restart/TCO95-CORE3/`. This is separate from the
  run's own restart directory and is easy to forget.
- The runscript itself,
  `~/esm_tools/runscripts/awiesm3/develop/awiesm3-develop-levante-TCO95L91-CORE3_PICAL_ccnice.yaml`.
  It is deliberately **not** tracked in git, so a clone will not bring it.

Paths beginning `/work/` are levante. Nothing under them exists on albedo.

## The software stack, which you should not improvise

`esm_tools`: branch `feat/awiesm3-v3.4-co2`, at `da82c6db9` when this was written. It carries
the v3.5 CORE3 defaults this run depends on: `cavity_gamma_scale: 0.6`, the sea-ice albedos
(`albsn` 0.80, `albsnm` 0.65, `albi` 0.70, `albim` 0.68, `albpnd` 0.20), `use_momix: false`,
`RCCNSEAICE: 15.0`, `RCL_INPPMIN: 70000.0`, and the mixing scheme with its IDEMIX forcing
paths. The runscript sets `install_esm_tools_branch` to this branch; if you change it the
run will not reproduce.

Components other than FESOM are fine as `esm_master` installs them. Each was verified on
2026-09-25 to sit exactly at its remote tip:

| component | branch |
|---|---|
| oifs-48r1 | `awiesm3-develop` |
| lpj_guess | `feat/lpjg-land-grid-from-coastline-4.1.11` (pins `guess.4111_lineA_13a4a06`) |
| xios | `main` |
| rnfmap | `enthalpy_on_ocean_side` |
| oasis3mct | `local_combined_fixes` (SMHI) |

**FESOM is the exception, and this is the part that matters.** `esm_master
install-awiesm3-develop` resolves `fesom-2.7-main` to branch `feat/cavity-gamma-scale`. That
branch is 53 commits behind `main` and **does not work for this run**. It has
`cavity_gamma_scale` and `use_atm_ice_tskin`, but it is missing both IDEMIX fixes: three
`MPI_DOUBLE_PRECISION` reductions remain in `gen_modules_cvmix_idemix.F90` and there is no
`nf90_fill_float` fallback in `gen_modules_read_NetCDF.F90`.

What you need is `main` + FESOM PR #1060 + FESOM PR #1087. Verified on 2026-09-25 to merge
with one trivial conflict, in `src/gen_modules_config.F90`, where `main` and #1060 each
appended a name to the same `namelist /run_config/` continuation line. Keep both
(`cavity_gamma_scale` and `use_atm_ice_tskin`).

Neither PR was merged when this was written. Check first: if both have landed, point
`configs/components/fesom/fesom-2.7.yaml`'s `2.7-main` at `main` and you are done.
`feat/cavity-gamma-scale` has outlived its purpose either way.

## Two things that will kill your first run

**1. `h0min` in the runscript.** The runscript carries

```yaml
    add_namelist_changes:
        namelist.ice:
            ice_therm:
                h0min: 0.5           # back to the original hard-coded value
```

0.5 is the code default, so this line does nothing at all. But `h0min` is only in the
`namelist /ice_therm/` group on the abandoned `feat/leadclose-h0min` branch. On `main`, on
`feat/cavity-gamma-scale` and on the integration branch it is a plain local variable in
`ice_thermo_cpl.F90`, not a namelist name, and FESOM aborts at the namelist read. **Delete
the line.** Do not try to keep it.

**2. The IDEMIX forcing files.** `awiesm3.yaml` now reads both from the FESOM pool:

```
${fesom.pool_dir}/forcing/idemix/fourier_smooth_2005_cfsr_inert_rgrid.nc
${fesom.pool_dir}/forcing/idemix/forcing_idemix_final_bin/FIN_STORMTIDE2_M2_plus_NYCANDER_CnoM2_bin_0.40deg.nc
```

with `fesom.pool_dir` = `${general.pool_dir}/fesom2/`, and `general.pool_dir` on albedo is
`/albedo/pool/`. Both files must exist under that layout before the first submit. The paths
the stock FESOM 2.7 namelist ships are useless: its surface file is a DKRZ `/pool/data/`
path and its bottom file is an `/albedo/pool/FESOM/fesom2.0/` path, so neither machine
resolves both.

## The gate: check this before you let a leg run

IDEMIX is only meaningful on a FESOM carrying #1087. Within the first minute of the compute
log you must see

```
--> IDEMIX total srf. energy Etot_srf =  0.3211645      TW
--> IDEMIX total tidal energy Etot_bot =  0.9135067     TW
--> IDEMIX Etot_bot after normalizing =  0.9179869      TW
```

0.9179869 is the forcing file's own `Etot_in_TW_bin` attribute, so this is a real check and
not a self-consistency one. If you see `Etot_srf = 6.05e-04` and `Etot_bot = -1.54e16`
instead, you have loaded a FESOM without the fixes. **Kill the job.** That exact signature
NaN'd the ocean in 100 seconds on 2026-09-24 and cost a 47-node slot.

## Facts you should not have to rediscover

- **The staged library can be a symlink, and `md5sum` follows it.** `<exp>/bin/fesom/lib/fesom/libfesom.so`
  is what the run loads, not the path in `bin_sources`. esm_tools keeps per-leg snapshots
  `libfesom.so_<leg>` and points `libfesom.so` at one of them. On 2026-09-24 it pointed at a
  pre-fix snapshot while every checksum looked right. Use `ls -l` before `md5sum`. The
  `fesom` executable is a driver stub and is byte-identical across all these builds, so a
  matching binary checksum proves nothing. Only the `.so` distinguishes them.
- **Chained legs do not read your runscript edits.** `<exp>/scripts/<name>.yaml` is a symlink
  into a per-leg archive. Edit the source and submit with `-U`, or the change is silently
  ignored.
- **The date file advances past a crashed leg.** `<exp>/scripts/*.date` is not evidence of
  progress. Count years in `outdata/` and restart files instead.
- **Preserve failed run directories.** Move them aside with a descriptive name. Never
  `rm -rf`. Deleting one has already destroyed the only evidence for a claim that was then
  challenged.
- **Never run two `esm_master` in one `model_codes` directory.** Concurrent installs
  overwrite each other's shared `comp-*_script.sh`.
- **OASIS `rmp_*.nc` weights belong to one FESOM partition.** Changing `fesom.nproc` is safe
  because the lookup is keyed on it, but an `oasis3mct.ini_restart_dir` carrying `rmp` files
  will shadow that lookup and smuggle wrong-partition weights in silently.
- **Contained-run mode:** `general.use_venv: True` is equivalent to `--contained-run`, but
  `--open-run` in `original_command` overrides the yaml and is self-perpetuating through the
  chain. Break the chain and resubmit fresh if you need to change it.
- Do not submit `esm_runscripts` from inside a Slurm job. Inherited `SLURM_MEM_PER_NODE`
  kills `srun`.

## Running it

`nyear: 10` and `restart_rate: 10` must match, or the leg cannot finish. On levante the run
used 47 nodes and reached about 77 SYPD under TKE+IDEMIX (3:18:32 for 2080-89, 2:58:12 for
2090-99, against about 86 SYPD for `cvmix_TKE` alone and 75 for KPP). Albedo's node layout
and per-node core count differ, so treat those numbers as a levante baseline and measure
your own before planning a campaign.

Evaluation lives in the campaign repo on levante, at
`/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive`:
`scripts/analysis/ccnice_regional_drift.py`, `ice_volume_series.py`,
`giomas_volume_targets.py`, `coupled_annual_store.py`, and `scripts/figures/`. The full
`reval` tool is `/work/ab0246/a270092/software/release_evaluation_tool2`, driven by a
per-window config under its `configs/`; copy `awiesm3_v35_pical_ccnice_2090s.py` and change
`model_version`, `pi_ctrl_start`/`end`, `spinup_end` and `out_path`.

## Things not to do

- **`surf_relax_s` stays 0.** Sea surface salinity restoring is excluded as a lever, even for
  Weddell convection. This is not negotiable and has been settled repeatedly.
- **No hemispheric tuning parameters.** `albsn` and `albsnm` are global. The reason one acts
  mainly in the south and the other mainly in the north is measured surface state
  (`d(alb)/d(albsn)` is 0.908 in SH DJF against 0.149 in NH JJA), not a hemispheric switch.
- **Price any albedo change in W/m2 first.** The Arctic runs away above about +5.5 W/m2; the
  Antarctic is stable to at least +10.4. The sensitivities are state-dependent and have
  changed by a factor of three as the Arctic cooled, so re-measure rather than reusing a
  number from a note.
- **Do not open pull requests or merge requests unprompted** in any component repository.
  Branches, commits and pushes are fine. Jan controls what goes back upstream.

## What to watch, and the question this spin-up has not answered

NH April volume is the number that decides whether `cvmix_TKE+cvmix_IDEMIX` keeps its place.
It was +1.145 +- 0.968 per decade over 2080-2099 and sits 3.37 above the early GIOMAS record.
If it keeps climbing, the scheme comes back out.

The upper-ocean heat drain is the mechanism behind the cold drift and nobody has tried to
stop it. It is worth knowing whether it is the mixing scheme, the Southern Ocean freshwater
regime, or both, before another lever is proposed.

The Southern Ocean deficit is not a vertical-mixing problem. Two regimes have been
characterised (salty/venting and fresh/warm CDW) and `spp`, `h0`, MLE, stronger GM and
basal-melt 0.6 were all tested in 2026-09 without fixing it in 20 years.
