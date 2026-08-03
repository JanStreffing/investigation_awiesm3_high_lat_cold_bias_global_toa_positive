# Two coupled tests: G4, and G4 + snow depletion

Both are namelist-only on one binary. Copy your `09C_06V_CRUNCEPinit_newSeaIce` runscript
twice and change nothing but the `oifs: add_namelist_changes:` block.

1. **Code — plain fast-forward**, your `f3ccacb` is the merge-base (5 ahead, 0 behind):
   `cd /work/ab0246/$USER/model_codes/awiesm3-develop/oifs-48r1` then
   `git fetch origin investigation/high-lat-cold-bias-round10 && git merge --ff-only origin/investigation/high-lat-cold-bias-round10`
2. All five commits are **no-ops at their defaults**, so the merge alone cannot change 09C.
3. **Rebuild** `esm_master recomp-awiesm3-develop/oifs`, and check the `libsurf.SP.so` md5
   actually changed before submitting.
4. **Test A (G4)** — keep your existing four namelist groups exactly as they are, add the
   two `NAMCLDP` lines and the whole `NAMSURFTUNE` group:
   ```yaml
           fort.4:
               NAMCLDP:
                   RVICE: 0.16                     # keep
                   RCL_INPSEA: 0.2                 # + D2b, Southern Ocean
                   RCL_INPPMIN: 70000.0            # + D2b, below ~700 hPa only
               NAMGWWMS:                           # keep as is
                   GGAUSSB: -0.5
               NAEPHY:                             # keep as is
                   LRDALB: False
               NAMCUMF:                            # keep as is
                   ENTSTPC3: 1
               NAMSURFTUNE:                        # + entire group is new
                   "ECE_TUNE_RVRSMIN(3)": 1000.0   # evergreen needleleaf
                   "ECE_TUNE_RVRSMIN(4)": 1000.0   # deciduous needleleaf
                   "ECE_TUNE_RVRSMIN(9)": 225.0    # tundra (80 -> 225)
   ```
5. **Test B (G4 + snow depletion)** — identical, plus two more lines in `NAMSURFTUNE`:
   ```yaml
                   ECE_SNOW_SCF: 1                 # sub-grid depletion on
                   ECE_SNOW_SCF_Z0: 0.016          # calibrated vs satellite
   ```
6. Quote the `RVRSMIN` keys — the brackets need them in YAML. `ECE_SNOW_SCF*` take no
   brackets and no quotes.
7. It must be **`NAMSURFTUNE`, not `NAMECECFG`**: `NAMECECFG` is read twice (arpifs
   `ecearth.F90` and surf `surfece.F90`) and an unknown name there aborts the arpifs read
   at `su0yoma.F90:152`.
8. **Verify it took** — `grep -a SURFECE_APPLY_TUNING NODE.001_01` must show the three
   `RVRSMIN` lines (A and B) and the `snow cover scheme 1 (Niu&Yang)` line (B only).
   Printed at setup, so a 1-day run is enough to check before committing 100 y.
9. `LRDALB: False` stays — the snow-cover fraction feeds the albedo blend through that
   branch, so B only does anything with it as you already have it.
10. Watch first: snow melt-out date and May/June snow cover, then tree cover. **B matters
    for the forest** — G4 alone delays melt-out a further 4 days, and later melt is later
    growing-season onset for LPJG.

AMIP counterparts: `amip_G4_tundra`, `amip_I1_scf` in `/work/bb1469/a270092/runtime/oifsamip-cy48/`.
Reasoning in `report/report.pdf` §"Route B resolved" and `notes/RUNS_AND_PARAMETERS.md`.
