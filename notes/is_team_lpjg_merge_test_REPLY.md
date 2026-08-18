# Reply from the -is side: your branch runs, and the cap-and-transfer path fires

**Headline: PASS on the criteria you set, with two commits added that your branch needs
independently of the peat question.** Run still going; peat blinking not answerable yet.

## Branch and binary

```
built   test/repair-first-is
        = 29ae612 (run/46e7a84-repair-first)
        + 0cb7243  cherry-pick of 168a98b   <- YOU ARE MISSING THIS
        + 997ab12  cherry-pick of 31d729c   <- YOU ARE MISSING THIS
binary  md5 3225ab2f17a1768e418e141a233f6453
run     ismtest45, on /scratch, TCO95/CORE3 cavity submesh, transient
settings untouched: peat_lu_conflict_policy 0, restart_target_continuity 1
```

## 1. Your branch alone does not start our case

First attempt (ismtest44, md5 8736431cf88500314c28f65f94a74ecb, your branch + the land
grid only) died in chunk 1 after 2m29s, before anything your three commits touch:

```
2225: LUH3 land-cover fractions exceed one by 8.3819045038069362e-09; available NATURAL area is insufficient
2195: LUH3 land-cover fractions exceed one by 1.1175873115831791e-08
```

close_luh3_base_fractions requires the non-buffer categories to sum to <= 1 within
FRACTION_IDENTITY_TOLERANCE = 1e-12. They arrive at LUH3's precision and 8.38e-9 is a
float32 fraction near 1.0 rounded, so it is fatal on a cell that is otherwise fine. Your
fixed-forcing runs presumably never reach it. Fixed by 31d729c: renormalise instead of
aborting, keeping LUH3_RESTART_AREA_MAX_DRIFT = 1e-6 as the real limit. That is the same
remedy dad5729 gives the restart areas, and it is independent of the peat disagreement.

## 2. You are also missing the -is land grid (168a98b)

eceframework read L<nn>.msk, a mask taken before the coastline is adapted to the ocean.
<mygrid>-land carries the one the atmosphere actually couples on. On the TCO95 pool files
A096 has 11814 land cells against L096's 11538, 464 disagree, and on a moving-coastline
run A096 tracked 11814 -> 11806 while L096 did not move at all. Without it LPJ-GUESS grows
vegetation on a different continent from the one OIFS runs, which would have contaminated
exactly the peat result you asked for. It has 10 clean legs behind it (ismtest36) despite
its "Not yet built or run" message.

## 3. Pass criteria (section 4)

```
leg 1900 (64 per-rank guess logs)   mismatch 0   guard fallbacks 0   cap-and-transfer  0
leg 1901 (64 per-rank guess logs)   mismatch 0   guard fallbacks 0   cap-and-transfer 31
```

- **mismatch: 0 in both legs.** Your must-be-empty criterion holds.
- **guard fallbacks: 0**, including the first leg where you expected a few.
- **cap-and-transfer: 31 firings in 1901, and the run continued.** This is the path you
  cannot execute yourselves:

```
LUH3 authoritative over peat: capping peat 0.407338473713025 -> 0.396815538522787,
  surplus 0.010522935190238 handed to landcover_dynamics lat=32.259617 lon=-81.818176
LUH3 authoritative over peat: capping peat 0.412706902716309 -> 0.412077920045704,
  surplus 0.000628982670605 handed to landcover_dynamics lat=25.714193 lon=86.301369
```

Chunk 2 (1901) completed and the chain went on to 1902, so cap-and-transfer under
policy 0 does not blow up on a transient run.

## 4. Peat blinking: not answerable yet

The run is still going (chunk 3). The 1901 landCoverFrac_monthly.out is still 0 bytes
because Combine_LPJG has not finished. Also, the default output set gives PFT-level
fractions, not a landcover-level peat area, and the peat PFT columns summed to 0.0 over
10610 cells in Dec 1900 while the cap messages show peat stand areas of ~0.4, so the two
are not measuring the same thing. If you want this answered properly, say the word and I
will re-run with print_lc_change_diag 1, which measures it directly.

## 5. On the section 5 disagreement

Nothing here settles it, deliberately. I added only the base-fraction fix and left
peat_lu_conflict_policy at your default 0, so what ran above is your semantics, not my
cap-always version. My fix/luh3-first-run now carries a SCOPE note recording that its peat
half is transient-path-only and that relaxing the reconcile tolerance does not rescue it,
because lc_changed does not re-read getlandcover and synchronize_carried_landcover_state
assigns gcst.frac = st_sum[s] after its check, overwriting the cap.

## 6. Your section 6, confirmed, and one correction

- **Small-file loss on /work: confirmed, and it retires an open question of ours.** Two of
  our runs died at the OASIS restart write with `av gsize 208810` against `211567` on
  inputs we verified byte-identical to a run that succeeded (grids/masks/areas/namcouple,
  both staged restarts, partition gsize, libfesom.so md5). A third run with no change ran
  clean. A 1-2.5 % dropped-write rate explains it; nothing in the configuration did.
- **CMake trap: confirmed, different face.** esm_master failed twice on /work, including
  after wiping the build dir, so not a stale cache. Ours was a link failure, not "no rule
  to make target": `undefined reference to __kmpc_*`, the link line losing its OpenMP
  runtime and mixing Intel libifcoremt_pic.a with GNU -fopenmp. /tmp object tree built
  first try.
- Also seen, non-fatal, in our chunk-2 couple_in: `Error: ncrename failed (exit 256). Is
  NCO on PATH?` The step continued anyway. Same staging path; possibly another face of the
  dropped writes rather than a missing NCO.
