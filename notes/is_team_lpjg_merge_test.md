# For the -is side: please test our LPJ-GUESS branch before we merge

**Ask:** build `run/46e7a84-repair-first` and run one `ismtest`-style transient case on it.
If it clears, we merge our branch and yours into `lpj_guess_awiesm3` together.
**Why you:** the branch changes the peat/LUH3 path, and your transient runs exercise a
code path our fixed-forcing runs cannot reach at all.

---

## 1. What to build

```
repo   git.smhi.se/e8891/lpjg-4.1
branch run/46e7a84-repair-first        head 29ae612
base   lpj_guess_awiesm3 (dad5729) + 3 commits -- a clean fast-forward, nothing diverged
```

Three commits on top of quasi-main:

| commit | what |
|---|---|
| `8c068f3` | repair-first bookkeeping checks; `peat_lu_conflict_policy` switch |
| `ab45fdd` | difference land-cover targets against the carried value, not a re-sum of the stands |
| `29ae612` | latch at stand-type granularity; match the check's summation order |

## 2. Settings — defaults are what we want tested

Change nothing. Both new switches default correctly for your case:

```
peat_lu_conflict_policy     0   (default) LUH3 authoritative -- CMIP7-correct
restart_target_continuity   1   (default) the restart fix
```

Do **not** set `peat_lu_conflict_policy 1`. That is for fixed-forcing control runs only.

If you want the land-cover diagnostics, add `print_lc_change_diag 1` — one line per
gridcell per year when anything changes. Verbose; optional.

## 3. What we specifically need from your side

Under `fixed_LU`, policy 0 **refuses** rather than capping peat, because the transfer
machinery does not run there — `landcover_init` returns before `landcover_dynamics`. So
we have never executed the cap-and-transfer path. **Your transient runs are the only
place it runs.** That is the gap your test closes.

Your `fix/luh3-first-run` reaches the same place by a different route (cap always, no
switch). We deliberately have not merged either into the other — see §5.

## 4. Pass criteria

Your existing one: **`ismtest39` gets past 1901.** Plus these greps over the per-rank
`guess*.log`:

```bash
grep -rh "LUH3 previous LC/ST mismatch\|LUH3 current LC/ST mismatch" <run>/work/run*/guess*.log
grep -rh "Carried stand-type metadata disagrees"                     <run>/work/run*/guess*.log
grep -rh "LUH3 authoritative over peat"                              <run>/work/run*/guess*.log
```

- **mismatch — must be empty.** Two of our three commits exist because it was not. If it
  reappears, stop and send us the line; it means a summation order we have not found.
- **guard fallbacks — expect a few on the FIRST leg only** (restart from a foreign parent
  configuration, which is what the guard is for). Fallbacks on a later self-restart are a
  finding, please report.
- **"LUH3 authoritative over peat"** — this is the cap-and-transfer firing. Seeing it, and
  the run continuing, is the positive result we cannot produce ourselves.

Also worth a look: does peat area stay sane across leg boundaries, or blink? At HR we
measured a 9.5 % peat-cell swing driven by a 0.33 % cropland wobble; that is the symptom
this branch is meant to remove.

## 5. The one real disagreement, not to be resolved by merging

Both branches rewrite the same block in `framework/externalinput.cpp`. Git will merge
them; the semantics will not.

- **Yours:** cap peat, always, and let `landcover_dynamics` move the surplus.
- **Ours:** a switch — LUH3-authoritative by default; refuse under `fixed_LU` because the
  transfer cannot run there; physical-peat-authoritative as an opt-in for control runs.

Ours refuses where yours caps, in exactly the configuration where capping silently
rewrites bookkeeping without moving the stand. We measured that on 2026-08-14:
`st=2 metadata=0.0886 physical=0.1917`, caught one call later by the next check.

Neither of us should merge over the other. Once your test passes we should agree one
combined behaviour and land it in one commit.

## 6. Two environment traps that will cost you a day

Both bit us hard on 2026-08-17/18.

- **`/work` drops small-file writes at ~1-2.5 %.** Measured: 2/200 on `/work/ab0246`,
  5/200 on `/work/bb1469`, **0/2000 on `/scratch`**. Large-file I/O (cmorization) never
  notices; a ~5000-file `esm_runscripts` staging copy essentially always fails, with
  `shutil.Error`/`FileNotFoundError` on files that exist and read fine. It also corrupts
  `guess.ins` during staging, which surfaces as a misleading **"Bad instruction file!"**
  and every LPJ-GUESS rank exiting 99 at startup.
  **Run on `/scratch`.** Ours submitted first try there after ~20 failures on `/work`.
- **CMake builds fail the same way on `/work`** (`link.txt` / `compiler_depend.ts` "no
  rule to make target"). **Build with the object tree on `/tmp`.**

## 7. Reporting back

Enough for us to act on:

- branch/commit built, and binary md5
- did `ismtest39` pass 1901
- the three grep counts from §4
- any peat blinking across leg boundaries

Background: `notes/lpjg_restart_design.md` (why the restart changes exist) and campaign
report §Round 29, both in
`JanStreffing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive`.
