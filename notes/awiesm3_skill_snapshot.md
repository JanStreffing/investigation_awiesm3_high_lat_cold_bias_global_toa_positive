---
name: awiesm3
description: Working on AWI-ESM3 / AWI-CM3 (OpenIFS 48r1 + FESOM2 + LPJ-GUESS + OASIS + XIOS) - where runtime truth is recorded, how to check a claim before spending a 48-year run, which observational references are independent, and the recurring failure modes that have cost this project the most. Use before proposing a tuning lever, diagnosing a crash, asserting what a parameter is, or comparing the model against an observation.
---

# Working on AWI-ESM3

This model is expensive to ask questions of. A 48-year AMIP run is ~9.5 h; a coupled
leg is longer. A grep is seconds. **The dominant cost in this project has not been
compute, it has been running simulations to answer questions that were already
answered somewhere on disk.**

Of 38 claims one campaign had to retract, roughly **45 % were answerable at a desk in
seconds and were not checked**, at a cost of about six wasted 48-year runs. Every one
of them came from resolving a claim at the wrong layer.

## 1. Resolve the claim at the layer where it is decided

This is the whole skill in one table. Answering anywhere else is guessing with extra
steps.

| the claim | the only place that settles it |
|---|---|
| what a parameter **is** | the run's `fort.4` / namelist — **never** the source default |
| **why** it has that value | `git log -S`, runscript comments, the project issues |
| whether a lever was **already tried** | the run inventory (`preflight.py`), not memory |
| how far a run **got** | `ifs.stat` — one line per timestep |
| why a run **died** | `NODE.001_01`, `drhook*`, and *every* rank — not just rank 1 |
| whether an observation is **independent** | what produced it (see §4) |
| whether a **null** means anything | the detection threshold, computed *first* |
| whether a **mechanism** works that way | the code path, traced end to end |
| how to **rebuild** after a source change | `esm_master recomp-` (§4), not the incremental `comp-` script |

### Source defaults are not runtime values

A default assignment in a setup routine (`sucldp.F90`, `susrad_mod.F90`, `su0phy.F90`…)
is what applies **when the namelist is silent**. It is not what the model ran with.

Two real failures from this pattern:

- `sucldp.F90:628` sets `NAERCLD=0 ! 0 = no aerosol interactions`, which was reported as
  "the aerosol–cloud path is off and CCN is a uniform 125". The runtime `fort.4` carries
  `LMACV2SP_CCNF = .true.` — a *second* CCN path, active, making `PCCN` a field.
- "`RVRSMIN` 500 is the IFS default and was halved in this build" — a misread of the
  table. `git log` showed both values present in the initial 48r1 commit.

**Always** confirm against `fort.4` before asserting a value.

## 2. Where runtime truth lives (the work folder)

Full reference: <https://awi-esm3.readthedocs.io/en/latest/workfolder.html>. The main
esm_tools log is **not** in `work/` — it is in the run's `log/` folder.

**OpenIFS**
- `fort.4` — every OpenIFS namelist, including `NAERAD` and `NAMORB`. **The authority
  on what the atmosphere ran with.**
- `NODE.001_01` — the detailed run log. Setup routines print resolved parameter values
  here (though some `WRITE`s are conditional, so absence is not proof).
- `ifs.stat` — one line per timestep. First thing to read for "how far did it get".
- `drhook*` — traceback when enabled. `gstats.xml`, `meminfo.txt` for timing/memory.

**FESOM2** — `namelist.config`, `.oce`, `.ice`, `.icepack`, `.io`, `.dyn`, `.tra`,
`.cvmix`. `fesom.clock` (first number = last written timestep) for restart state.

**OASIS** — `namcouple` is the coupling contract: fields, remapping, frequencies.
`debug.root.{rank}` / `nout.{rank}` for what it actually parsed and exchanged.
`grids.nc`, `masks.nc`, `areas.nc` when remapping looks wrong — a conservative remap
that silently returns zero is usually a mask problem.

**XIOS** — `file_def.xml` decides what output exists at all. Check it before concluding
a variable is missing. `xios_client_{rank}.out` for write failures.

**LPJ-GUESS** — `guess.log`, `*.ins` (`run_coupled_*.ins` is the entry point),
`LPJ-GUESS_monthlyoutput.txt` for what vegetation actually handed to OpenIFS.

### Crash triage order

1. `ifs.stat` — where it stopped.
2. **Grep every rank**, not rank 1. Search for `ABOR1`, `MPL_ABORT`, `SIGSEGV`,
   `forrtl`, and the model's own error strings. A campaign once invented an
   infrastructure story because rank 1 was searched for tokens that do not match how
   OpenIFS aborts; the real abort was on `PROC=48`.
3. Only then consider the batch system.

## 3. Before proposing a lever: run `preflight.py`

`preflight.py PARAM` (in this skill's `scripts/`, and in the campaign repo under `scripts/analysis/`) reports the source default, namelist
reachability, **every run that has ever set it and to what value**, runscripts
(including withdrawn), `git log -S` provenance, and prior mentions in the report.

It exists because of three specific failures:

- **A duplicate.** `RCL_OVERLAPLIQICE 0.65→0.35` was submitted as new. It is `A1b`,
  already run 44 years; `A1a` is the same knob at `0.1`. The label `ovl` meant the
  liquid/ice **deposition** overlap and was read as radiative cloud overlap.
  **A label is not an identifier — the namelist value is.**
- **A stale note.** "The land-albedo residual is untouched by all 41 runs" survived in
  an open list after `K1` had already closed 74 % of it.
- **An un-asked why.** Lowering `RVICE` was proposed to warm the tropics, without
  checking that it had been *raised* from 0.13 to 0.16 precisely to suppress Southern
  Ocean high cloud.

If a parameter is not in a namelist, exposing it means declaration + association +
namelist entry. In OpenIFS `sucldp.F90` the namelist-settable variables use a module
`POINTER` bound in the "Associate pointers for variables in namelist" block; variables
bound only in the `ASSOCIATE` block are **default-setting only** and cannot be reached.
Move, do not duplicate — a name cannot be both.

**`NAMSURFTUNE`, not `NAMECECFG`.** `NAMECECFG` is read twice from the same `fort.4`,
and a Fortran namelist read dies with *"invalid reference to variable"* on any name the
reading module does not declare. That hard failure is also a useful *positive* test: if
a run starts and passes initialisation, the namelist accepted every name in it.

## 4. Rebuilding after a source change: use `esm_master recomp-`

```
cd /work/<proj>/<user>/model_codes
esm_master recomp-<setup>-<version>        # e.g. recomp-oifsamip-cy48
```

`esm_master` with no arguments lists every setup and the operations it supports
(`comp clean get update status log install recomp`).

**`recomp-` is `conf-` + `clean-` + `comp-`** — a *clean* rebuild. The per-setup
`comp-<model>_script.sh` is the *incremental* build. Use `recomp-` when a module's
declarations change, or when a failure is not understood; use `comp-` when the cause is
known and confined.

### `error #6404: This name does not have a type` — check accessibility first

The OpenIFS modules declare `PRIVATE` as the default and then export an explicit list:

```fortran
    PRIVATE
    ...
    PUBLIC ECE_CLIMR
    PUBLIC ECE_CLIMR_DMS      ! <- a new module variable is INVISIBLE without this
```

A new module-level variable therefore compiles fine and **still fails at every `USE`
site** until it is added to the `PUBLIC` list. This looks exactly like a stale-module
or build-ordering problem and is not one — a full `recomp-` will not fix it.

**Diagnose it with `nm`, not by guessing at the build system.** If the symbol is present
in the module's object file, the source and the compile are fine and the problem is
accessibility:

```
nm build/.../ecearth.F90.o | grep -i my_new_var
# ecearth_mp_ece_dms_ccn_sens_   <- present, so it is a PUBLIC problem, not a build one
```

(Real case, 2026-08-09: three rebuilds were spent on a stale-`.mod` theory — including a
full clean `recomp-` — for a variable that was simply never exported. `nm` would have
answered it in seconds.)

### Verifying a rebuild actually reached the run

The executable often does **not** change, because the physics lives in
`libarpifs.SP.so` and is linked dynamically — check the *library* timestamp, not the
binary. And esm_runscripts **stages a per-run copy** of the library at submit time, so
a rebuild cannot disturb jobs already queued or running: check
`<run>/work/lib/oifs/libarpifs.SP.so` to see what a given leg actually loaded. This is
also why a mid-campaign rebuild is safe, and why a *new* run picks up the change while
a running multi-leg job does not.

`strings` on the library is **not** a reliable check that a namelist name was exposed:
names appear there for other reasons (`RAMID` and `RCLCRIT` both show up without being
in `NAMCLDP`). The dependable test is functional — a Fortran namelist read aborts on an
undeclared name, so a run that starts and passes initialisation has accepted every name
in its `fort.4`.

## 5. Observational references: which are independent

The single most repeated error class is treating a non-independent reference as
independent.

- **ERA5 is IFS.** For an IFS-derived model it is a sibling, not an observation. It
  reproduced the model's Southern Ocean cloud bias almost exactly (low area *and* low
  reflection, in proportion). Use it for dynamics and thermodynamics; **never** as the
  arbiter on cloud, snow, soil or turbulent flux for this model.
- **Reanalysis-derived forcing is not truth.** CRUNCEP3 was ~21 W/m² brighter than
  satellite, which manufactured a "33 W/m² Siberian SW deficit" that did not exist.
- **CERES** is a genuine independent radiometer. But `cldarea` is a MODIS **mask**
  (counts optically thin cloud) and is *not* comparable to a model `tcc` (a radiative
  cover). It is reliable over dark ocean and unreliable over ice — CERES and MODIS
  disagree by 7.7 pp at 90S–65S.
- **Compare what the cloud DOES, not how much is counted.** CRE is a difference of two
  broadband fluxes, identically defined across model, reanalysis and satellite. Cloud
  area is not one quantity.
- **Period-match, or measure the offset.** A PI run against present-day CERES/ERA5
  carries an epoch offset that has been large enough to invert a conclusion (a tropical
  TOA "error" of −2.51 became −0.67 period-clean). Prefer the model's own present-day
  arm.

## 6. Before believing a result

- **State the detection threshold first.** Compute it from interannual scatter, then
  compare. A lever once read `+0.502 K` on a 4-year window and `−0.038 K` over 44 —
  pure noise, after two runs had been built on it. A 10-year coupled pair resolves only
  ~±0.75 K of Siberian DJF soil, so a null there is not evidence of safety.
- **Check the control is actually a control.** Two "land-only" levers were used as a
  null for a Southern Ocean diagnostic; both *contained* the ocean lever being tested.
- **Score on the quantity you are claiming.** CRE cannot distinguish cloud amount from
  cloud opacity. Mean snow cover cannot see a distributional change that mean-preserves.
  TOA cannot see a surface parameterisation error worth 20 K of soil temperature.
- **Averaging does not commute with a nonlinear operator.** `SCF(mean d) ≠ mean(SCF(d))`
  voided one analysis outright. Apply the function per sample, then average.
- **IFS TOA and surface fluxes are accumulated J/m² over the output step.** Divide by
  the accumulation period or every number is ~3600× too large.
- **Snow enthalpy** (`sf × 3.3355e8`) belongs in the surface budget and is worth
  ~1 W/m² — omitting it manufactures a spurious TOA-vs-surface gap.

## 7. Screen at one year before spending 44

A 1-year AMIP leg is ~11 minutes; the 44-year campaign standard is ~9.5 hours. Project
issue #170 found the cy48 high-cloud offset already present in year 1 while T2m and the
low-cloud dip lagged, and concluded 1-year runs suffice for testing tuning parameters.
**Measured on Southern Ocean SW CRE, that holds — with a threshold.**

Control interannual scatter is 0.712 W/m², so a pair of single years resolves
**±1.97 W/m²** against ±0.30 for a 44-year pair. Levers above that reproduce their
44-year answer at one year (A1a −6.56 vs −6.44; D2a −5.23 vs −4.80; D2b −2.01 vs −2.61).
Levers below it do not, and **B3 comes out with the opposite sign**.

So: **screen the radiative response at one year, never adopt on one.** Temperature is
explicitly excluded — the seasonal thresholds (DJF ±0.588 K) are 44-year numbers, and
#170's own lead/lag result is that T2m lags the cloud change. Compute the 1-year
threshold for your metric from the control's interannual scatter before trusting a short
run, exactly as for the long ones.

## 8. Check the prior record before designing a round

Beyond this repository: the **AWI-ESM/project_management issues** carry years of tuning
history that is not in any code comment — why a parameter has a non-default value, which
levers were tried and rejected, and which biases are known to be structural. Two of the
errors above were already answered there.

Search the issues, the campaign report, and `notes/` before proposing anything. If a
lever looks obviously good and nobody has tried it, the most likely explanation is that
somebody has.
