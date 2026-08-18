# Review of `lpjg_restart_design.md`

**Reviewer:** Claude Opus 5. **Date:** 2026-08-17.
**Source verified against:** `lpj_guess_repairfirst` @ `8c068f3` (working tree), all citations re-read.
**Verdict:** the diagnosis is real and §1–§3 hold on inspection. §4 — the section the argument
rests on — is false. That error propagates into the §6 measurement design and into Option A's
guard. Two cheaper options are missing from §7.

---

## 1. What was verified and holds

| Claim | Status |
|---|---|
| §2 `Landcover::serialize` is `arch & frac` only | confirmed, `guess.cpp:2966` |
| §2 `Gridcellst::serialize` carries `frac`, `frac_old_orig`; not `frac_old`/`frac_change` | confirmed, `guess.cpp:2902` |
| §2 `Stand::serialize` carries `frac`, pfts, patches, `first_year` | confirmed, `guess.cpp:1578` |
| §2 `Gridcell::serialize` carries `landcover` and every `st[i]` | confirmed, `guess.cpp:3521,3556` |
| §3 re-derivation block, verbatim including the renormalisation | confirmed, `externalinput.cpp:790-880` |
| §3 `lc.frac` zeroed and re-read from LUH3 as the new target | confirmed |
| §4 point 1: `getlandcover` runs twice per year under `fixed_LU` | confirmed, `framework.cpp:4481` |
| §7B: a layout change is read as garbage, not rejected | confirmed, `archive.h:27-39,115-124` |

The `archive.h` finding is the strongest thing in the document and it is exactly right.
`ArchiveStream` is `transfer(char*, std::streamsize)`; `operator&` writes
`sizeof(data)` raw bytes with no tag. There is no magic number, no version int, no field
name, and no length check on anything but `std::vector`. Adding a field to any serialiser
silently reinterprets every existing restart file from that offset onward. That argument
needs no further support and should survive any rewrite of the rest.

---

## 2. Blocking: §4 is false

`landcover_reconcile_restart_state` (`landcover.cpp:119`) **is**
`synchronize_carried_landcover_state` — it forwards to it at `landcover.cpp:129`:

```cpp
bool landcover_reconcile_restart_state(Gridcell& gridcell) {
	if(gridcell.nbr_stands() == 0) return false;
	return synchronize_carried_landcover_state(gridcell);   // :129
}
```

So the function §4 marks **"NOT REACHED"** at `4298` is reached every year at
`landcover.cpp:161`, on the restart path, from inside `landcover_init`. The §4 diagram
already prints `landcover.cpp:161  reconcile; return` and then draws the opposite
conclusion from it.

And the synchronisation is not conditional on the `>tol` branch. `landcover.cpp:3603`
and `:3625` assign unconditionally, with physical stands as truth:

```cpp
gcst.frac = st_sum[s];
gcst.frac_change = gcst.frac - gcst.frac_old;
...
gridcell.landcover.frac[lc] = lc_sum[lc];
gridcell.landcover.frac_change[lc] = ... - ...frac_old[lc];
```

The `>tol` branches above them only decide whether a `dprintf` fires. So the invariant
`frac == sum(stands)` **is** re-established, every year, in the direction "stands win".

§4's headline — *"why the invariant is not re-established"* — and its consequence 2
— *"the reconciliation is skipped exactly when nothing changes"* — are both the
opposite of what the code does. Consequence 1 (double `getlandcover`) is correct and is
the more important structural fact; it deserves to be the section's subject.

---

## 3. §6's premise is answerable statically, and is probably false

The untested premise is:

> that the serialised `frac` equals `sum(stands)` at write time, for a *self*-restart.

Tracing one `fixed_LU` year end to end settles the mechanism without a run:

1. `framework.cpp:4492` → `landcover_init` (every year, because `needs_landcover_init`
   includes `ECEARTH && fixedLUafter >= 0`)
   - `landcover.cpp:156` `getlandcover` **#1**: `frac_old ← physical`, `frac ← LUH3 target`
   - `landcover.cpp:161` `synchronize`: `frac ← st_sum` (**physical**), `frac_change ← 0`
2. `framework.cpp:718` → `landcover_dynamics` → `lc_changed`
   - `landcover.cpp:3366` `getlandcover` **#2**: `frac_old ← physical`, `frac ← LUH3 target`
   - `no_changes` → return at `3851` / `3927` / `3938`, before `4298`
3. End of leg, serialise: `frac` = **LUH3 target**, `sum(stands)` = **physical**

These are two different quantities — a target and a state — that coincide only after a
completed transition. Expect them to differ. §6's table then mis-reads that outcome:

> | differ structurally | the re-derivation is load-bearing; Option A is wrong |

That inference does not follow. They can differ structurally *and* Option A still be the
right fix, because A's value is not physical consistency — it is **target continuity**.
Under `fixed_LU` the LUH3 target is identical year to year, so `frac_old ← deserialised
frac` yields `frac_change == 0` **exactly**, `no_changes` holds, and nothing flips. The
0.33 % cropland wobble goes away because the comparison becomes target-to-target instead
of target-to-resummed-physical.

The measurement as specified cannot separate "differs because the parent configuration
differs" (reason 2) from "differs because `frac` is a target and `sum(stands)` is a state"
(structural, present on every self-restart). It will return the same verdict for both.

**Recommended replacement for §6.** Do not test the premise; it is settled. Measure the
*magnitude* instead, which is the part that governs the guard threshold:

- log `lc.frac`, `gcst.frac`, `sum(stands)` at **four** points per year — deserialise,
  after `getlandcover` #1, after `synchronize`, after `getlandcover` #2 — otherwise the
  drift cannot be attributed to a stage;
- branch the dprintf on `arch.save()` so serialise and deserialise are distinguishable
  (they are the same function);
- run **≥3 legs**, not 2. Leg 1 reads the 2000-year offline spin-up — the foreign parent
  of reason 2. Only the leg-2→leg-3 boundary is a clean self-restart. A 2-leg run gives
  one boundary and it is the contaminated one.

---

## 4. Option A does not survive as drafted

Two separate defects, both fixable, neither mentioned.

**4.1 The guard threshold is self-defeating.** A proposes: *"validate against `sum(stands)`;
on disagreement beyond float epsilon, fall back to re-derivation with a loud message."*
Per §3 above, disagreement beyond float epsilon is the **normal** state of every cell on
every restart. The guard fires everywhere, the fallback runs everywhere, and A degenerates
into D with extra logging. The guard must be a structural sanity bound —
`LC_METADATA_REPAIR_LIMIT` scale, the same band `synchronize` already treats as
"structurally sane but worth shouting about" — not float epsilon.

The §8 note inherits this: *"under A the right comparison becomes exact equality, and any
deviation is a bug rather than noise to threshold."* Exact equality against
`sum(stands)` is unreachable. Exact equality against **last year's target** is both
reachable and the thing A actually guarantees. That is the invariant worth asserting.

**4.2 The latch.** `getlandcover` runs twice per year, and `synchronize` overwrites
`frac` with `st_sum` in between (`landcover.cpp:3603`). By call #2 — the call that feeds
`lc_changed`, the one that decides `no_changes`, the only one that matters — the
deserialised `frac` is gone from memory. An implementation that reads `lc.frac` when it
happens to run is a **silent no-op on the load-bearing call**. A must latch the value at
deserialise time and consume it once, then mark it spent.

Neither defect is fatal to A. Both are fatal to A *as written*, and 4.2 in particular
would produce a change that compiles, runs, logs nothing unusual, and does not fix the
wobble.

---

## 5. Two options missing from §7

**E — stop calling `landcover_init` every year.**

```cpp
bool needs_landcover_init = (gridcell.nbr_stands() == 0) ||
                            (ECEARTH && ecearth.fixedLUafter >= 0);   // framework.cpp:4481
```

The comment immediately above it (`framework.cpp:4483-4488`) insists that a loaded
transient state *"must be left physically untouched here"* — and then the `fixedLUafter`
disjunct routes every `fixed_LU` restart into `landcover_init` → `synchronize` →
`stand.set_gridcell_fraction(... * scale)` at `landcover.cpp:3561`, which rescales every
stand in the cell. The code contradicts its own comment.

Dropping the `fixedLUafter` disjunct removes re-derivation #1 and the read-path stand
rescale outright, leaving the `landcover_dynamics` path — which is where the `-is` team's
`fix/luh3-first-run` already expects the work to happen (§9). This is plausibly smaller
than A and removes a cause rather than reinterpreting one. §4 assembled all the evidence
for it and §7 does not list it. It should be priced against A before either is written.

**B′ — sidecar file.** Write `frac_old`/`frac_change` to a small companion file beside
the archive rather than into it. Additive; backward compatible by construction (a missing
sidecar means an old restart, which falls back to re-derivation — exactly reason 2's
case); needs no magic number, no version int, and no v0 reader. It obtains B's correctness
without C's serialiser surgery and without invalidating a single existing file. C remains
the textbook answer, but on this fork's actual constraints B′ dominates it.

---

## 6. Smaller corrections

**§5's "complete list" is not complete.** It misses `framework/guess.cpp:882`
(`Stand::init_stand_lu`). That is stand *creation*, not movement of an existing
fraction, so §5's conclusion survives — but a section whose entire force is exhaustiveness
cannot have an unacknowledged gap. State the exclusion rather than omitting the line.

**§5's conclusion has a third case.** *"They can only disagree if they already disagreed
when written, or if the read path introduces the difference."* Line `3561` is inside
`synchronize`, which runs at read (`:161`). So the read path does not merely introduce a
difference in the metadata — it **moves the stands themselves**. That is a materially
different statement and it strengthens rather than weakens the overall case.

**`frac_old_orig` is serialised and never discussed.** §2 lists it in the table and §2's
prose then says only the deltas are missing. It is set from `gcst.frac` at
`externalinput.cpp:1645` and `:1806` and restored at `:1749`. Either it is a usable
checkpoint of the previous target — in which case A is nearly free and needs no new field
— or it is not, and §2 should say why in one line. Leaving it in the table unexamined is
the one place the document looks like it stopped early.

**§4's line numbers for the early returns** (`3851`, `3927`, `3938`) are right to within
a line or two in the current tree; worth re-pinning after any edit, since they are the
only citations in the document that are not exact.

---

## 7. Suggested disposition

Keep §1, §2, §3, §7B, §7C and §9 as they stand.

- Rewrite §4 around the double `getlandcover` call, and correct the reconcile claim: the
  invariant *is* re-established at `:161`, with stands as truth.
- Replace §6: drop the premise test, measure the per-leg drift magnitude at four points
  over ≥3 legs, and use it to set the guard threshold rather than to choose between options.
- Re-found Option A on target continuity, widen its guard to `LC_METADATA_REPAIR_LIMIT`
  scale, and specify the latch.
- Add E and B′ to §7. Price E against A first — E may be the smallest true fix, and unlike
  A it does not depend on what the archive happens to contain.
- Fix §5's list, and resolve `frac_old_orig` in one line.

The recommendation in §8 — do the measurement, then A, and keep C separate from A — is
sound as a *shape* even though both of its inputs need repair. Nothing found here argues
for D.

---

## 8. Note on method

This review was done by hand; there is no software-design-review skill in this
environment (`code-review`, `simplify` and `security-review` are all diff- and
bug-oriented and none of them reads a design document). The checks that actually found
things were mechanical enough to encode, should one be written:

1. Re-read every `file:line` citation rather than trusting it.
2. Grep every claim of the form "complete list" / "nothing else does X".
3. Follow the call graph for every "NOT REACHED" — here, the unreached function was
   reached under a different name three lines from a call site the document itself printed.
4. Trace one full cycle of the state machine and ask what each quantity *means* at each
   point, not just what it equals.
5. Test whether the decision table's outcomes actually discriminate its hypotheses.
