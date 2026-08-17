# LPJ-GUESS restart: the land-cover state is checkpointed as a target, not a state

**Status:** static analysis complete. The premise v1 wanted to measure is settled on
inspection; what remains to measure is a magnitude, not a yes/no (§6).
**Date:** 2026-08-17 (rev 2). **Source:** `lpj_guess_repairfirst` @ `8c068f3` + working tree.
**Rev 2 corrects a false claim in rev 1** — see §10.

---

## 1. The problem in one paragraph

A restart should be a checkpoint. This one checkpoints the pool-bearing stands
faithfully, and checkpoints the land-cover bookkeeping in a form that cannot be used as
a checkpoint: what lands in the archive is the **LUH3 target** from the last
`getlandcover` call of the year, while the stands hold the **physical state**. On read,
the target is discarded and `frac_old` is rebuilt by re-summing the stands. That
summation carries float roundoff, and the difference then meets equality-and-threshold
tests that convert it into whole stand types switching on and off. Measured at HR: a
0.33 % wobble in cropland area across leg boundaries, driving a 9.5 % swing in
peat-bearing cells (43 467 ↔ 48 049), with no carbon signature.

---

## 2. What is persisted

| Object | Serialised fields | File |
|---|---|---|
| `Stand` | `first_year`, **`frac`**, pfts, patches | `guess.cpp:1578` |
| `Gridcellst` | **`frac`**, `frac_old_orig`, `nstands`, `distinterval_st`, `diam_cut_low`, `nfert` | `guess.cpp:2902` |
| `Landcover` | **`frac`** | `guess.cpp:2966` |
| `Gridcell` | `landcover`, and `st[i]` for all `i` | `guess.cpp:3521,3556` |

**Not serialised:** `lc.frac_old`, `lc.frac_change`, `gcst.frac_old`, `gcst.frac_change`.

**`frac_old_orig` is not a usable previous-target checkpoint.** All three uses are
CROPLAND-scoped scratch inside the crop-fraction machinery — saved before a rescale
(`externalinput.cpp:1645`, `:1806`), restored as a fallback when a year's crop data is
missing (`:1749`). It does not cover other land-cover types and is not maintained as a
year-boundary quantity.

---

## 3. What is re-derived on read

`externalinput.cpp:790-880`, in `getlandcover`, on the `!is_cold_start` branch:

```cpp
// frac_old is not serialized, and serialized LC/ST metadata can come from a
// parent configuration with different enabled land covers. Reconstruct the
// previous state from physical stands before reading the new target.
    physical_st[stand.stid]      += area;      // re-sum every stand
    physical_lc[stand.landcover] += area;
    physical_total               += area;
...
    lc.frac[i]     = physical_lc[i];
    lc.frac_old[i] = physical_lc[i];           // <-- the load-bearing assignment
    gcst.frac      = physical_st[i];
    gcst.frac_old  = physical_st[i];
```

`lc.frac` is then zeroed and re-read from LUH3 as the new target. The code names its own
cost a few lines down: *"physical_total is a sum over every stand in the cell, so it
carries whatever double [roundoff]"*.

Both stated reasons are real. Reason 1 (`frac_old` not serialised) is true. Reason 2
(parent configuration may differ) is genuine and is exactly the leg-1-from-spin-up case
— but it is paid on **every** restart rather than handled as the exception.

---

## 4. Control flow: two re-derivations, and what `frac` means at serialise

```
guess_coupled          framework.cpp:4492   landcover_init(...)         <-- pass 1
  landcover_init       landcover.cpp:156      getlandcover() #1           frac_old <- physical
                                                                          frac     <- LUH3 target
                       landcover.cpp:161      reconcile  ==  synchronize  frac     <- st_sum
                       landcover.cpp:163      return
simulate_day_ece       framework.cpp:718    landcover_dynamics(...)     <-- pass 2
  landcover_dynamics   landcover.cpp:3844     lc_changed()
    lc_changed         landcover.cpp:3366       getlandcover() #2         frac_old <- physical
                                                                          frac     <- LUH3 target
                       landcover.cpp:3851     if(no_changes) return
```

**`getlandcover` runs twice per year** under `fixed_LU`, because `needs_landcover_init`
is `(nbr_stands()==0) || (ECEARTH && fixedLUafter >= 0)` (`framework.cpp:4481`). The
re-derivation therefore happens twice, and the *second* one is the load-bearing call —
it is what feeds `lc_changed` and decides `no_changes`.

**The invariant is re-established, in the direction "stands win."**
`landcover_reconcile_restart_state` **is** `synchronize_carried_landcover_state`; it
forwards at `landcover.cpp:129`. Its assignments are unconditional —
`gcst.frac = st_sum[s]` (`:3603`), `landcover.frac[lc] = lc_sum[lc]` (`:3625`) — with the
`>tol` branches above them gating only a `dprintf`.

**Consequence — what is actually in the archive, and it depends on the path.** The last
write to `lc.frac` decides it:

| path | last write to `frac` | archive holds |
|---|---|---|
| `fixed_LU`, any year | `getlandcover` #2 | **LUH3 target** |
| transient, year with change | `synchronize` @ `4298` | **achieved physical** |
| transient, no-change year | `getlandcover` #2 | LUH3 target |

`landcover_init` does not run at all in a transient run: `needs_landcover_init` is false
once `fixedLUafter < 0` and stands exist (`framework.cpp:4481`, and the comment above it
says so explicitly). So transient years take pass 2 only, and a year with real change
runs to `4298`, where `frac ← st_sum`.

Target and state coincide only just after a completed transition. Expecting them to be
equal in general is a category error, and it is the one rev 1 made.

---

## 5. Who may move a stand fraction

`set_gridcell_fraction` callers:

```
modules/landcover.cpp: 429, 500, 511, 571, 771, 835     land-cover change machinery
modules/landcover.cpp: 3561                              our renormalisation (new, rev-2 tree)
framework/guess.cpp:   882                               Stand::init_stand_lu -- CREATION
```

Excluding creation, every writer is inside land-cover change. Two things follow:

1. `frac` and the stands cannot drift apart *during* a year through some third path.
2. The read path does not merely rewrite metadata — line `3561` is inside `synchronize`,
   which runs at read (`:161`), so on our tree **the read path rescales the stands
   themselves.** That strengthens the case rather than weakening it.

---

## 6. What to measure

The v1 premise — "does serialised `frac` equal `sum(stands)`?" — is answered in §4:
**no, structurally, on every restart**, because one is a target and the other a state.
No run is needed for that, and a run designed to test it would return the same verdict
for the structural case and for the reason-2 parent-configuration case, which are
different problems.

What is *not* known is the **magnitude** of the target-to-physical gap on a clean
self-restart, which is what sets the guard threshold in §7A.

Instrumented run:

- log `lc.frac`, `gcst.frac`, `sum(stands)` at **four points per year** — at
  deserialise, after `getlandcover` #1, after `synchronize`, after `getlandcover` #2 —
  otherwise drift cannot be attributed to a stage;
- branch the dprintf on `arch.save()`, since serialise and deserialise are the same
  function;
- **≥3 legs.** Leg 1 reads the 2000-year offline spin-up, which is reason 2's foreign
  parent; only the leg-2→leg-3 boundary is a clean self-restart. A 2-leg run samples one
  boundary and it is the contaminated one.
- **count `force_small_change` activations per year** (`landcover.cpp:3374`). If A drives
  this to zero without stranding real area, that is a second argument for A rather than a
  caveat against it — see §7A.

Existing evidence for scale, all from contaminated leg-1 restarts: stand-sum drift
3.3e-9 to 4.2e-7; metadata-vs-physical up to 3.3e-2 at 45.35 N 53.65 E.

---

## 7. Options

**A — Latch the checkpointed target and use it as `frac_old`.**
Value is **target continuity**, not physical consistency. Under `fixed_LU` the LUH3
target is identical year to year, so `frac_old ← last year's target` gives
`frac_change == 0` *exactly*, `no_changes` holds, and nothing flips. The wobble
disappears because the comparison becomes target-to-target instead of
target-to-resummed-physical. Requirements:

- **Latch at deserialise.** By `getlandcover` #2 the deserialised value has been
  overwritten twice (by #1 and by `synchronize`). Capture it on read, consume it once,
  mark it spent. An implementation that reads `lc.frac` wherever it happens to run is a
  **silent no-op on the load-bearing call** — it compiles, runs, logs nothing, and does
  not fix the wobble.
- **Structural guard, not float epsilon.** Target and physical differ beyond epsilon on
  every cell of every restart (§4), so an epsilon guard fires everywhere and A degenerates
  into D with extra logging. Use a `LC_METADATA_REPAIR_LIMIT`-scale bound — the band
  `synchronize` already treats as structurally sane but worth shouting about.
- **Transient runs: mostly a non-issue, and possibly an improvement.** The worry was that
  differencing target-to-target would strand any gap between intended and achieved area.
  Per §4 it largely does not arise: in a transient year *with* change the archive already
  holds the achieved physical state, so A reproduces today's semantics bit-exactly rather
  than by re-summation. The exposure is one narrow case — a cell whose transition
  under-delivered, followed by a year with no change there.

  And in that case today's behaviour is not obviously the desirable one. Under-delivery
  is deliberate: `NEW_ST_THRESHOLD` and `filtered_amount` drop sub-resolution areas on
  purpose. A persistent target/physical gap at a near-zero category then trips

  ```cpp
  if(fabs(gcst.frac_change) > eps && (gcst.frac <= eps || gcst.frac_old <= eps))
      force_small_change = true;                       // landcover.cpp:3374
  ```

  which overrides the no-change early return (`:3426`) and runs the full land-cover
  machinery every year to re-attempt a transfer that will be filtered again for the same
  reason. That is not a correction channel; it is churn, and plausibly a contributor to
  the wobble. A stops it.

  **Genuinely open, and nobody appears to have decided it:** is a deliberately-filtered
  sub-resolution area meant to be retried forever, or written off? The code currently
  does the former by accident.

**E — Stop calling `landcover_init` every year.**
Drop the `fixedLUafter` disjunct from `needs_landcover_init` (`framework.cpp:4481`). The
comment directly above it insists a loaded transient state *"must be left physically
untouched here"*, and the disjunct then routes every `fixed_LU` restart into
`landcover_init` → `synchronize` → the stand rescale at `landcover.cpp:3561`. The code
contradicts its own comment. Removes re-derivation #1 and the read-path stand rescale
outright, leaving the `landcover_dynamics` path — where the `-is` team's
`fix/luh3-first-run` already expects the work to happen. Plausibly smaller than A, and it
removes a cause rather than reinterpreting one. **Price this against A first.**

**B′ — Sidecar file.**
Write `frac_old`/`frac_change` to a small companion file beside the archive. Additive and
backward compatible by construction: a missing sidecar means an old restart, which falls
back to re-derivation — reason 2's case exactly. No magic number, no version int, no v0
reader. Obtains B's correctness without C's surgery and without invalidating a single
existing file.

**B — Serialise the deltas into the archive.**
Strictly more correct, and **breaks every existing restart file silently**.
`ArchiveStream` is a raw byte stream (`transfer(char*, streamsize)`, `archive.h:38`);
`operator&` writes `sizeof(data)` raw bytes with no tag, no magic, no version, no field
name. A layout change is read as garbage rather than rejected. Do not do this without C.

**C — Version the archive, then B.**
Magic + version int per element, readers branch. The textbook answer. On this fork's
constraints B′ dominates it: same correctness, none of the serialiser surgery.

**D — Status quo plus `LC_METADATA_REPAIR_LIMIT`.**
Works today. Restart-continuation stays non-reproducible and every future check in this
area inherits the same noise floor.

---

## 8. Recommendation

1. **Price E against A.** E may be the smallest true fix, and unlike A it does not depend
   on what the archive happens to contain or on getting a latch right. It is a one-line
   change to a disjunct whose own adjacent comment argues against it.
2. If E is rejected, do **A** with the latch and the structural guard, and settle the
   transient-run question first.
3. Run the §6 instrumentation regardless — it sets the guard threshold for A and gives
   the before/after evidence for E.
4. **B′ over C** if the deltas are ever wanted. Keep either well away from A.

Nothing found here argues for D.

Note for whoever picks this up: `lc_bookkeeping_tol() = INPUT_RESOLUTION * 1e-3 = 1e-9`
(`landcover.cpp:80`) is scaled off *input* precision but applied to a sum of
restart-state floats. Under A the assertable invariant is exact equality against **last
year's target** — reachable — not against `sum(stands)`, which is not.

---

## 9. Related

- Peat-vs-LUH3 conflict and the two-policy switch: `run/46e7a84-repair-first`, `8c068f3`.
- Campaign report §Round 29 (`report/report.tex`).
- The `-is` team's `fix/luh3-first-run` caps peat and relies on the transfer in
  `landcover_dynamics`. Correct in the transient path; under `fixed_LU` the
  `landcover_init` reconcile at `landcover.cpp:161` fires first — measured 2026-08-14 as
  `st=2 metadata=0.0886 physical=0.1917`.

---

## 10. Corrections to rev 1

Kept for honesty; rev 1 is in git history at `8322a84`.

- **§4 claimed the end-of-year reconciliation is "NOT REACHED" in no-change years, and
  that the invariant is therefore never re-established. False.**
  `landcover_reconcile_restart_state` forwards to `synchronize_carried_landcover_state`
  at `landcover.cpp:129`, so it runs every year at `:161`. Rev 1's own diagram printed
  that line and drew the opposite conclusion. Caught in review.
- **§6's premise test was unanswerable as designed.** Target and state differ
  structurally on every restart, so the test could not separate that from the
  parent-configuration case. Replaced with a magnitude measurement.
- **§7A's float-epsilon guard was self-defeating** for the same reason, and A lacked the
  latch, without which it is a silent no-op.
- **§5's "complete list" omitted `Stand::init_stand_lu` (`guess.cpp:882`).** Creation
  rather than movement, so the conclusion stood, but the omission was not stated.
- **`frac_old_orig` was tabled and never examined.** Now resolved: CROPLAND-scoped
  scratch, not a previous-target checkpoint, so A is not free.
- **Options E and B′ were missing**, and E is a candidate for the smallest true fix.
