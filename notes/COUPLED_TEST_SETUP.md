# Two coupled tests: G4, and G4 + snow depletion

Everything is namelist-driven, so both tests share one binary. Base your runscript on
`Tuning_test_09C_06V_CRUNCEPinit_newSeaIce`, change only the namelist blocks below.

1. **Get the code.** In `/work/ab0246/$USER/model_codes/awiesm3-develop/oifs-48r1`
   (branch `movcav-landice+co2-concdriven`, currently at `f3ccacb`), cherry-pick four
   commits from `investigation/high-lat-cold-bias-round10` — none of them are in your
   branch yet, and they are all no-ops at their defaults:
   ```bash
   git fetch origin investigation/high-lat-cold-bias-round10
   git cherry-pick 2630bd1 1004cba 28b5542 9f12d79
   ```
   `2630bd1` adds `RCL_INPSEA`/`RCL_INPPMIN` (needed for D2b), `1004cba` adds the
   `&NAMSURFTUNE` namelist, `28b5542` adds the snow-depletion scheme, `9f12d79` sets its
   calibrated default. Expect a conflict only if you have touched `cloudsc.F90`.

2. **Rebuild** with `esm_master recomp-awiesm3-develop/oifs` and **check the library md5
   actually changed** before submitting — a stale object has cost us five runs before.

3. **Test A — G4** (the best AMIP configuration: boreal +0.952 K, SO SW RMSE 6.88→4.80):
   ```yaml
   oifs:
       add_namelist_changes:
           fort.4:
               NAMCLDP:
                   RCL_INPSEA: 0.2
                   RCL_INPPMIN: 70000.0
               NAMSURFTUNE:
                   "ECE_TUNE_RVRSMIN(3)": 1000.0
                   "ECE_TUNE_RVRSMIN(4)": 1000.0
                   "ECE_TUNE_RVRSMIN(9)": 225.0
   ```

4. **Test B — G4 + snow depletion.** Identical, plus two lines in the same
   `NAMSURFTUNE` block:
   ```yaml
                   ECE_SNOW_SCF: 1
                   ECE_SNOW_SCF_Z0: 0.016
   ```

5. **`NAMSURFTUNE`, not `NAMECECFG`.** `NAMECECFG` is read twice (arpifs `ecearth.F90`
   and surf `surfece.F90`); putting these entries there aborts the arpifs read at
   `su0yoma.F90:152`. Keep `RVICE: 0.16` in `NAMCLDP` if your baseline had it.

6. **Check it took.** `NODE.001_01` must show, once per run:
   ```
   SURFECE_APPLY_TUNING: RVRSMIN(  3)   250.0000 ->   1000.0000
   SURFECE_APPLY_TUNING: snow cover scheme  1 (Niu&Yang) z0= 0.0160 ...
   ```
   If those lines are absent the namelist never reached the model. A 1-day run is enough
   to check — the message is printed at setup, before the first timestep.

7. **Why B might matter more than A for the forest.** G4 buys AMIP summer temperature but
   *delays* snow melt-out by a further 4 days (24→28 May against 11 May observed). Later
   melt is later growing-season onset for LPJ-GUESS, so G4 alone could improve T2m while
   still starving establishment. Test B removes that: it fixes the depletion so melt-out
   lands within a day of satellite.

8. **What to look at first**, beyond the usual: snow-cover melt-out date and May/June
   snow cover, then whether tree cover establishes. Both tests should be run long enough
   for vegetation to respond; the AMIP result cannot tell us whether 48 % of the bias is
   enough.

9. **Expect B to differ from A mainly in spring.** The scheme is inert in midwinter and
   does nothing for the autumn cover bias (Sep/Oct stay too snowy — a separate problem).

10. **AMIP counterparts for comparison** are `amip_G4_tundra` and `amip_I1_scf` under
    `/work/bb1469/a270092/runtime/oifsamip-cy48/`. Full reasoning is in
    `report/report.pdf` §"Route B resolved" and `notes/RUNS_AND_PARAMETERS.md`.
