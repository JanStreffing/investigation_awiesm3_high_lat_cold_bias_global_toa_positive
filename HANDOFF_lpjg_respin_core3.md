# Handoff: rebuild the LPJ-GUESS forcing and spin-up for the new CORE3 mesh

Written 2026-08-26. Starting point for whoever picks this up; not a complete recipe.
Everything below was established today and is worth not re-deriving.

## Why this is needed

The 40-year mesh test `11X` (new CORE3, otherwise identical to `11W`) crash-loops. The
failure is **not** the mesh and not the atmosphere. LPJ-GUESS ranks exit 99 on:

```
Fixed peat invariant failed in after coupled restart load at lon 54.782608, lat 72.466919,
year 1350: map=0.0086111099999999999 physical=0 LC=0 ST=0; totals physical=1 LC=1 ST=1.
State was not repaired.
```

That check is new in LPJ-GUESS **tag 4.1.11** (`8832661`, Laszlo Hajdu, 2026-08-24), whose
commit message promises to "validate peat, land fractions, stand structure, and soil state
after restart loading". The peat map says the cell is 0.86 % peat, the restart carries none,
4.1.11 attempts repair, cannot close it, and calls `fail()`.

It is gated only by `if(!run_landcover || !run[PEATLAND]) return;`
(`framework/externalinput.cpp:4371`). Both are on, so there is no namelist switch, and
relaxing the check is the wrong move: this campaign has already decided once that widening a
tolerance does not restore comparability.

**So the blocker is the STATE, not the code.** The state in use is
`LPJG-SPINUP_2000Y_TCO95_CORE3_CRUNCEPandCERES_daily_variability/.../lpjg_state_3901`,
produced long before 4.1.11 by a binary with no fixed-peat handling. The old binary
`guess.raupach_junebase_gated` reads it fine, which is how `11W` ran 50 years.

The new mesh compounds it: the coastline moves, so LPJ-GUESS also meets gridcells the
2000-year spin-up never simulated.

## The job

Two steps, in order.

**1. New 10-year AMIP run** with the latest atmospheric tuning, configured with the OIFS
output-control files that write exactly the fields LPJ-GUESS needs as forcing.

**2. New 2000-year LPJ-GUESS spin-up** driven by that forcing, using the LPJ-GUESS land-sea
mask that matches the new CORE3. That mask should already exist from ocp-tool; check before
generating one.

The product is a `lpjg_state_*` that 4.1.11 accepts and that matches the new coastline.
`11X` then reruns unchanged.

## Facts you should not have to rediscover

**Mesh discipline.** `core3` = 220509 nodes, `core3_beta` = 211567. The shipped mesh was
replaced in place on 2026-08-26 and the old one renamed; the fix concerned cavity depths.
`restart/CORE3` no longer exists. Analyse pre-2026-08-26 runs with `core3_beta`, later runs
with `core3`, and never mix: numpy will area-weight a 211567-node field against the first
211567 entries of a 220509-node array and return a plausible wrong number. Use
`scripts/analysis/meshguard.py`, which raises instead. The new `core3` ships **without**
`fesom.mesh.diag.nc`; reval needs that file and will have to have one generated.

**Binary.** Built and staged today, never successfully run:
`model_codes/awiesm3-develop/bin/guess.4111_landgrid_raupach`,
md5 `4b46d3564ab169ee52ccdbabfc5e485b`. From branch
`feat/lpjg-land-grid-from-coastline-4.1.11` = tag 4.1.11 plus the four commits unique to the
coastline branch (the `-land` gridlist, cold-start of a cell with no state, and the two
Raupach `GUE_Z0HV` commits). The three LUH3 tolerance workarounds and the peat
import/revert pair were dropped because 4.1.11 fixes the cause. **The branch is local only,
unpushed.** Build with `esm_master recomp-lpj_guess-4.1.2`, never cmake directly, and verify
the md5 changed.

**Latest atmospheric tuning** is LX4 + RSBLB with S4 removed. Runscripts under
`~/esm_tools/runscripts/oifsamip/`, closest arm `..._SB2_sblb2.yaml`. S4
(`RCL_INPPMIN` 50000) was retired on 2026-08-24: coupled it reverses sign on the metric it
was adopted for (Siberian JJA +0.590 K on removal, resolved) and `11V` gives the best
campaign CMPI at 0.7385. Do not reintroduce it.

**Leg length must be self-consistent.** esm_tools sizes the coupled leg from `nyear`, but
FESOM takes `run_length` from `restart_rate`. `nyear: 1` against `restart_rate: 10` builds a
one-year window while telling FESOM to run ten years; FESOM then reaches the end with no
restart written and the leg cannot continue. Keep `nyear`, `restart_rate` and
`restart_first` equal.

**A crashed leg can report success.** LPJ-GUESS ranks exiting 99 still left the job
`COMPLETED 0:0`, so esm_tools advanced the date file and resubmitted. `11X` produced three
legs in six minutes and zero years before it was caught. Check years-of-output, not exit
codes.

**Copies can truncate silently.** Staging `dist_1792` into `core3/` stopped at 159 of 3585
files with no error and no process running. `/work` is at 99 %. Verify file counts after any
large copy.

## Open question worth resolving first

Ask Laszlo whether 4.1.11 expects a regenerated state or carries a migration path for
pre-4.1.11 states. If a migration exists, both steps above become unnecessary. A cheap
discriminator meanwhile: run the 4.1.11 binary on the **old** mesh with the existing state.
If the invariant still fires, it is purely state compatibility and has nothing to do with
the cavity-depth question `11X` was built to answer.

## What `11X` is for, once unblocked

Isolating whether corrected cavity depths fix the Weddell polynya diagnosed on 2026-08-26:
Oct-Dec `MLD2` at 90-60S is 758.6 m against WOA18's 148.8 m, a quarter of the band convects
past 1000 m and the deepest node reaches 4160 m, on a surface 0.20 psu too salty. That
convection mines heat from 700-2000 m, where the model is 0.358 K too cold against PHC3,
melts ice, drops albedo and raises absorbed shortwave. It is the mechanism behind the
positive TOA and its drift. It will settle if MLD2 falls toward 148.8, the >1000 m area
fraction falls from 25 %, and the upper-50 m salinity bias falls from +0.181 psu.
