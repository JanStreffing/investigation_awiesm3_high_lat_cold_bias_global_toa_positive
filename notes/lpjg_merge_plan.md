# Merge plan: getting our LPJ-GUESS and esm_tools work back onto the main lines

**Date:** 2026-08-18. **Status:** nothing merged yet; conditions below are unmet.

Four divergent lines exist. They are not equally ready, and two of them conflict
semantically with each other. This plan says what merges, in what order, and what has
to be true first. **A condition that is not met blocks the merge — it is not advisory.**

---

## Current topology

`lpj_guess_awiesm3` (`dad5729`) is quasi-main.

| branch | ahead | contents | verdict |
|---|---|---|---|
| `feat/lpjg-land-grid-from-coastline` | +4 | gridlist maskfix `168a98b` | **merge first** |
| `run/46e7a84-repair-first` | +3 | repair-first checks, peat policy, restart continuity | **merge second, conditions below** |
| `origin/fix/luh3-first-run` (-is) | +3 | their peat cap, base-fraction renormalisation | **merge WITH ours, jointly** |
| `run/50yr-arms-raupach` | +3 | Raupach roughness | **do not merge** |

`run/46e7a84-repair-first` is a clean fast-forward: quasi-main has zero commits it lacks.
That makes it mechanically safe and says nothing about whether it is correct.

---

## Merge 1 — the gridlist maskfix (`168a98b`)

Build the LPJ-GUESS gridlist from `<mygrid>-land` rather than `L<nn>`.

**Why first:** best-evidenced thing we have. It carried 11I and 11J through 50 and 40
model years respectively. It is also currently stranded on a feature branch while the
production arms that validated it were built from it, which is the wrong way round.

**Conditions**

1. 11I complete at 50 years on a binary containing it. **MET** — 11I 1350-1399.
2. No LPJ-GUESS aborts attributable to it across those runs. **MET.**
3. Cold-start path present, since the new gridlist adds ~307 cells with no spin-up state
   and every one of them aborted without it. **MET** — `0218bbb`, already on quasi-main.
4. Someone other than us has built it at least once. **NOT MET** — covered by the -is
   test (`notes/is_team_lpjg_merge_test.md`), whose branch descends from quasi-main.

---

## Merge 2 — repair-first checks, peat policy, restart continuity

Three commits: `8c068f3`, `ab45fdd`, `29ae612`.

**Conditions**

1. **The A/B measurement exists.** Two runs, one binary, one namelist line apart:
   `restart_target_continuity` 1 vs 0, under `fixed_LU`, with `print_lc_change_diag 1`.
   Required result: **LCDIAG == 0 with the switch on, > 0 with it off.**
   *Status: UNMET. Attempt 1 (S2_fix/S2_ctl, 2026-08-18) produced no measurement at all --
   both arms reported LCDIAG 0, but every fpc.out and cpool.out was zero bytes, so neither
   arm simulated a year and the counter never reached the code that increments it. Cause was
   two more roundoff tolerances in externalinput.cpp that the landcover.cpp relaxation left
   behind: LUH3_RESTART_AREA_MAX_DRIFT at 1e-6 against an observed 9.28e-6, and Laszlo's
   close_luh3_base_fractions at 1e-12 against an observed 4.37e-8. Both now 1e-4 in c226bbf.
   Attempt 2 (S3_ctl / S4_fix, binary guess.restart_ab_roundoff md5 441c3f04) reached
   three legs but ALSO cannot produce the measurement: the LCDIAG instrument itself was
   blind. It sat ~40 lines above where change_st is accumulated and before
   change_gross_lcc exists, so it reported change_st as 0.0 unconditionally and could
   never fire on it -- and change_st is the stand-type term, which is exactly the
   granularity the defect lives at. Fixed in lpjg 4d91c3b by moving it to immediately
   before the no-change decision and reporting all four terms that decision tests.
   Binary guess.restart_ab_lcdiagfix md5 6022479c. Attempt 3 must use it.

   Attempt 2 is still worth finishing for conditions 4 and 5 (no LC/ST mismatch, guard
   fallbacks confined to leg 1) -- those do not depend on LCDIAG.

   Failed arms preserved: S2_*.failed_roundoff_20260818, and S3_fix (0-byte 44.state).*

   *THREE instrument failures preceded any data on this measurement: Slurm job state
   (esm_tools reports COMPLETED on a crashed leg), the log path (dprintf lands in
   work/**/guess*.log, not log/*.log), and variable ordering (above). Each produced a
   plausible zero. Before scoring attempt 3, verify the control arm produces LCDIAG > 0
   -- if it does not, suspect the instrument before concluding the arms agree.*

   *Do not score this from Slurm job state. esm_tools reports COMPLETED on a crashed leg and
   advances the date file -- S2_ctl walked forward two legs on restarts that were never
   written. Score on non-zero annual output.*
   If both arms come out 0 the test proved nothing and must be redesigned — that is a
   failure, not a pass.
2. **`peat_lu_conflict_policy 1` runs to completion at least once.** It has never
   executed: it was silently broken by a Jinja trailing-newline bug until 2026-08-17, and
   every run since has died before reaching it. *UNMET.*
3. **`peat_lu_conflict_policy 0`'s cap-and-transfer path runs at least once.** Impossible
   for us — under `fixed_LU` it refuses by design. Requires the -is transient test.
   *UNMET, delegated.*
4. **No `LUH3 previous/current LC/ST mismatch` in any test run.** Two of the three commits
   exist because this fired; it is the sharpest regression signal available.
   *Status: absent from the last scratch runs. Provisionally met, reconfirm on S2.*
5. **Guard fallbacks confined to the first leg.** `Carried stand-type metadata disagrees`
   is expected when restarting from a foreign parent configuration and nowhere else.
   A fallback on a self-restart means `LC_PREV_FRAC_MAX_DRIFT = 0.1` is mis-tuned.
   *UNMET — needs a >=2-leg run.*
6. **Agreement with -is on the peat conflict.** See below. *UNMET.*

Conditions 1, 2 and 5 are all satisfied by one 3-leg fixed_LU run per arm. Condition 3 is
the -is test. Condition 6 is a conversation.

---

## The peat conflict — resolve before either side merges

Both branches rewrite the same block of `framework/externalinput.cpp`.

- **-is:** cap peat to available NATURAL, always, and let `landcover_dynamics` move the
  surplus with its C/N/water.
- **ours:** a run-time switch. LUH3-authoritative by default (CMIP7 requires the land-use
  dataset be reproduced); **refuse** under `fixed_LU` because the transfer machinery does
  not run there; physical-peat-authoritative as an opt-in for control runs.

Git will merge these. The semantics will not. Ours refuses exactly where theirs caps, and
we have the measurement for why: capping under `fixed_LU` rewrites the bookkeeping while
the physical stand keeps its area, and the next check catches it —
`st=2 metadata=0.0886 physical=0.1917`, 2026-08-14.

Their reasoning is right for the transient path and ours is right for the fixed one.
The likely landing point is theirs plus our `fixed_LU` guard, as one commit, not two
merges. **Neither side merges over the other.**

---

## Do not merge: Raupach

`run/50yr-arms-raupach` (`1385d0c`). Falsified — the winter mechanism needs a canopy the
model does not grow: measured lambda_stem 0.0092 at 2.89 m, giving z0 = 0.0132 m against
IFS's 0.058 m, so it makes the winter surface 4.4x *smoother*. Keep the branch as a
record, gated off behind `ifraupachz0`. Campaign report §Round 29.

---

## esm_tools

```
awiesm3-is-orography-and-smb-coupling  vs  feat/awiesm3-v3.4-co2
  ahead 69   behind 5   merge-base 30469198f
```

**Do not merge wholesale.** The branch is mixed: ice-sheet/orography/SMB coupling plus the
AWI-ESM3 runscript and namelist work. A single merge drags the ice-sheet coupling into the
CO2 line.

**The set is far smaller than 69.** Classified 2026-08-18 by the paths each commit
touches, then checked by hand — most of the branch is iceberg / PISM / ice2fesom work.
The AWI-ESM3-relevant commits are **three**:

| commit | touches | note |
|---|---|---|
| `fe2c3f86f` | the `peatauth` template + 4 AWI-ESM3 runscripts | carries the coupling-field placement fix: `add_veg_atm_fields` in `lpj_guess`, `add_atm_veg_fields` in `oifs`, **not** `oasis3mct` — misplaced, `esm_runscripts` exits without creating the experiment directory |
| `c0755b6ea` | `configs/components/lpj_guess/lpj_guess.cmip.yaml` | fixed 1850 defaults |
| `8141777ac` | `namelists/lpj_guess/ecearth.ins.j2` | `state_remap` to the nearest gridcell when the land moves |

**Skip `1c3d60767`** (fesom xios instantaneous stress / daily m_ice): verified
byte-identical to upstream `540dc3518`, which landed independently. Cherry-picking it
would duplicate or conflict.

**Also note two of the five upstream commits we are behind interact with our campaign
directly:** `1bcb76563` "Split the ICMGG and LPJ-GUESS suffixes, and hold LPJ-GUESS on
the old map" and `ae28c3d5f` "Take the repaired ICMGG and slt on every version but the
released tags". The target branch pins LPJ-GUESS to the old soil map; 11I overrides that
with `lpjg_slt_suffix: "_v2"`. Whoever merges must decide whether the pin stays, given
11I measured the repair as a trade rather than free — Siberian DJF soil +0.308 K but the
ERA5 JJA bias worsening from -0.673 to -1.123. That is a science decision, not a merge
mechanic.

**Conditions**

1. Rebase or merge the 5 upstream commits in first. *UNMET.*
2. Cherry-pick list agreed and reviewed — a handful, not 69. *UNMET.*
3. No ice-sheet-only commit in the set. Check every hunk.
4. One runscript from the set actually submits and runs after the cherry-pick. *UNMET.*

---

## Environment preconditions for any of this

Both cost a day on 2026-08-17/18 and will cost the next person the same.

- **`/work` drops 1-2.5 % of small-file writes.** Measured 2/200 (`ab0246`), 5/200
  (`bb1469`), 0/2000 small files on `/scratch`. Coarse I/O never notices; a ~5000-file
  `esm_runscripts` staging copy essentially always fails. It also corrupts staged
  `guess.ins`, which surfaces as a misleading **"Bad instruction file!"** with every
  LPJ-GUESS rank exiting 99. **Run on `/scratch`.**
- **CMake fails on `/work`** (`link.txt` / `compiler_depend.ts` missing). **Build with the
  object tree on `/tmp`** — `BUILDDIR=/tmp/... comp_lpjg_worktree.sh`.

- **`/scratch` is ~100x cleaner but NOT clean.** Measured 2026-08-18 across every run
  staged there: **1 zero-byte `.state` in 9472 files, 0.011 %**. The one failure was a
  54 MB LPJ-GUESS restart file whose two source copies were byte-identical at
  54,518,971 bytes -- so it is a dropped write on the staging copy, not a truncated
  source, and it is not confined to small files. At ~64 state files per leg boundary
  that is ~0.7 % per boundary: rare, but a long multi-leg run will hit it. It killed
  S3_fix at leg 2 with "failed to read index for state file".

  **Defensive measure for any multi-leg run: verify the staged `.state` files are
  non-empty before the leg starts.** The failure otherwise surfaces one leg later as a
  model abort, and esm_tools crash-forwards past it.

Worth a DKRZ ticket: 200-file write-then-stat reproducer on /work, plus the /scratch
large-file case above (1/9472, with both source copies intact and md5-identical).

---

## Order

```
1. finish our A/B            (conditions 2.1, 2.2, 2.5)
2. -is transient test        (conditions 1.4, 2.3)
3. agree the peat behaviour  (condition 2.6)
4. merge 168a98b
5. merge repair-first + the agreed peat commit, together
6. esm_tools cherry-picks, separately
```

Raupach never merges. Everything above waits on step 1, which is one 3-leg run per arm.
