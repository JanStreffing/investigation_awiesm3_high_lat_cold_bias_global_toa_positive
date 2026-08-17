# LPJ-GUESS restart: the land-cover state is checkpointed and then discarded

**Status:** static analysis complete; one premise still needs a measurement (§6).
**Date:** 2026-08-17. **Source read:** `lpj_guess_repairfirst` @ `dad5729` + local changes.

---

## 1. The problem in one paragraph

A restart should be a checkpoint. This one is a checkpoint for the pool-bearing stands
and a *re-initialisation* for the land-cover bookkeeping laid over them. The bookkeeping
is written into the archive, read back, and then overwritten by a re-summation of the
stands. That re-summation is not bit-identical to the value it replaces, and the
difference meets equality-and-threshold tests a few lines later, which convert it into
whole stand types switching on and off. Measured consequence at HR: a 0.33 % wobble in
cropland area across leg boundaries, driving a 9.5 % swing in peat-bearing cells
(43 467 ↔ 48 049), with no carbon signature.

---

## 2. What is persisted

Everything below survives the round trip. The archive is complete for the stands and
*almost* complete for the bookkeeping.

| Object | Serialised fields | File |
|---|---|---|
| `Stand` | `first_year`, **`frac`**, pfts, patches | `guess.cpp:1578` |
| `Gridcellst` | **`frac`**, `frac_old_orig`, `nstands`, `distinterval_st`, `diam_cut_low`, `nfert` | `guess.cpp:2902` |
| `Landcover` | **`frac`** | `guess.cpp:2966` |
| `Gridcell` | `landcover`, and `st[i]` for all `i` | `guess.cpp:3521,3556` |

**Not serialised:** `lc.frac_old`, `lc.frac_change`, `gcst.frac_old`, `gcst.frac_change`.

So the *level* (`frac`) is checkpointed at all three levels. Only the *deltas* are not.

---

## 3. What is re-derived on read

`externalinput.cpp:793-870`, in `getlandcover`, on the `!is_cold_start` branch:

```cpp
// frac_old is not serialized, and serialized LC/ST metadata can come from a
// parent configuration with different enabled land covers. Reconstruct the
// previous state from physical stands before reading the new target.
    physical_st[stand.stid]      += area;      // re-sum every stand
    physical_lc[stand.landcover] += area;
    physical_total               += area;
...
    lc.frac[i]     = physical_lc[i];           // overwrite the restored value
    lc.frac_old[i] = physical_lc[i];
    gcst.frac      = physical_st[i];
    gcst.frac_old  = physical_st[i];
```

`lc.frac` is then immediately zeroed and re-read from LUH3 as the new *target*, so the
load-bearing assignment is `frac_old ← re-summed stands`. The comment gives two reasons.
Reason 1 (`frac_old` is not serialised) is true but does not require discarding `frac`,
which *is*. Reason 2 (parent configuration may differ) is a real concern, but it is
being paid on every ordinary restart rather than handled as the exception it is.

The code names the cost itself a few lines down: *"physical_total is a sum over every
stand in the cell, so it carries whatever double [roundoff]"*.

---

## 4. Control flow — why the invariant is not re-established

`frac == sum(stands)` would be forced at `landcover.cpp:4298`
(`synchronize_carried_landcover_state`, "a completed operation must leave physical stand
sums equal to the targets"). It is not reached in the ordinary case.

```
guess_coupled          framework.cpp:4492   landcover_init(...)          <-- first
  landcover_init       landcover.cpp:156      getlandcover()               re-derive #1
                       landcover.cpp:161      reconcile; return            <-- returns
simulate_day_ece       framework.cpp:718    landcover_dynamics(...)      <-- second
  landcover_dynamics   landcover.cpp:3844     lc_changed()
    lc_changed         landcover.cpp:3366       getlandcover()             re-derive #2
                       landcover.cpp:3849     if(no_changes ...) return    <-- returns
                       landcover.cpp:4298     synchronize_...              <-- NOT REACHED
```

Three separate no-change early returns (`3851`, `3927`, `3938`), all guarded on
`no_changes`, all before `4298`.

Two consequences:

1. **`getlandcover` runs twice per year** under `fixed_LU`, because
   `needs_landcover_init` is `(nbr_stands()==0) || (ECEARTH && fixedLUafter >= 0)`
   (`framework.cpp:4481`). The re-derivation happens twice, not once.
2. **The reconciliation is skipped exactly when nothing changes**, which is the normal
   year in a fixed-forcing control run. `frac` and `sum(stands)` are only forced equal
   in years that had a real land-cover change.

---

## 5. Who may move a stand fraction

Complete list of `set_gridcell_fraction` callers outside the accessor itself:

```
modules/landcover.cpp: 429, 500, 511, 571, 771, 835     land-cover change machinery
modules/landcover.cpp: 3561                              our renormalisation (new)
```

All inside land-cover change. **Nothing else in the model moves a stand fraction.** This
matters: `frac` and the stands cannot drift apart *during* a year. They can only
disagree if they already disagreed when written, or if the read path introduces the
difference — which it does.

---

## 6. What is established, and the one thing that is not

**Established statically:**

- `frac` is in the archive at all three levels and is discarded on read.
- The re-derivation is a float summation over stands; its error floor is float epsilon,
  not `INPUT_RESOLUTION`. Measured 3.3e-9 to 4.2e-7 against a 1e-9 tolerance.
- The reconciliation that would restore the invariant is skipped in no-change years.
- Only land-cover change moves stand fractions.
- Downstream, `gcst.frac_old == 0.0` with `new_frac < NEW_ST_THRESHOLD`
  (`externalinput.cpp:1090`) and `raw_luh3_crop_sum > lc.frac[CROPLAND]`
  (`crop_fraction_scale`) turn that difference into a discrete category flip.

**NOT established — the premise the cheap fix rests on:**

> that the serialised `frac` equals `sum(stands)` at write time, for a *self*-restart.

All measurements to date are leg-1 restarts from the 2000-year offline spin-up — a
foreign parent configuration, exactly reason 2 above — where metadata and stands
genuinely differ by up to 3.3e-2 (`st=0 metadata=0.4605 physical=0.4932`, 45.35 N
53.65 E). That says nothing about leg 2+.

**Measurement that settles it:** a 2-year, 2-leg TCO95 run with a dprintf of `lc.frac`,
`gcst.frac` and `sum(stands)` at serialise and at deserialise.

| outcome | meaning |
|---|---|
| identical | Option A is correct and trivial |
| differ at ~1e-16 | Option A works with an exact-equality-plus-epsilon guard |
| differ structurally | the re-derivation is load-bearing; Option A is wrong |

---

## 7. Options

**A — Trust the checkpoint (cheap, no format change).**
Use the deserialised `frac` as `frac_old` instead of re-summing. Validate against
`sum(stands)`; on disagreement beyond float epsilon, fall back to re-derivation with a
loud message — which is then exactly the parent-configuration case reason 2 was written
for. ~20 lines, no archive change, no invalidated restart files, revertible.
*Blocked on the §6 measurement.*

**B — Serialise the deltas too.**
Add `frac_old`/`frac_change` to `Landcover::serialize` and `Gridcellst::serialize`.
Strictly more correct. **Breaks every existing restart file**, including the spin-up
state, and does so *silently*: `ArchiveStream` is a raw byte stream
(`transfer(char*, streamsize)`, `archive.h:38`) with no magic number, no version tag and
no field names, so a layout change is read as garbage rather than rejected.

**C — Version the archive, then B.**
Magic + version int per element; readers branch. The right long-term design, and what
makes B safe. Cost: touches the serialiser core, needs a v0 reader for existing files,
and we would be the only testers.

**D — Leave it; keep thresholding the symptom.**
Status quo plus `LC_METADATA_REPAIR_LIMIT`. Works, but restart-continuation stays
non-reproducible and every future check in this area inherits the same noise floor.

---

## 8. Recommendation

Do the §6 measurement, then **A** if it supports it. It removes a cause rather than
adding a seventh workaround, and it is small enough to offer upstream rather than carry.

**C** is the correct end state and is more attractive than it first looks: the fork has
not merged upstream in ~18 months, so "every existing restart file" effectively means
*our* restart files, which we control. It should not be bundled with A — a silent
failure mode deserves its own change and its own testing.

Note for whoever picks this up: the `INPUT_RESOLUTION`-derived tolerance
(`lc_bookkeeping_tol() = INPUT_RESOLUTION * 1e-3 = 1e-9`, `landcover.cpp:80`) is only
defensible while the value is treated as freshly derived from input. Under A the right
comparison becomes exact equality, and any deviation is a bug rather than noise to
threshold.

---

## 9. Related

- Peat-vs-LUH3 conflict and the two-policy switch: branch `run/46e7a84-repair-first`,
  commit `8c068f3`.
- Repair-first bookkeeping limits: same commit, `modules/landcover.cpp`.
- Campaign report §Round 29 (`report/report.tex`).
- The `-is` team's `fix/luh3-first-run` caps peat and relies on the transfer at
  `landcover_dynamics`. That works in the transient path; under `fixed_LU` the
  `landcover_init` return at `landcover.cpp:161` aborts first — measured 2026-08-14 as
  `st=2 metadata=0.0886 physical=0.1917`.
