# REcoM in the OpenIFS-coupled build: surface pCO2 falls because sea-level pressure is zero

State of 2026-10-05, evening. The cause is established from the output and the source; the fix is drafted and built, and its one-month test is running. Nothing is pushed.

## Symptom

First seen in PICAL_cc_v350val (carbon-cycle variant, concentration-driven, 2140-2141), the first run that wrote REcoM diagnostics. PICAL_cc has identical ocean fields, so it has the same problem; it just wrote no BGC output.

Global-mean surface values, area weighted:

| month | xCO2atm (ppm) | pCO2s (µatm) | dpCO2s (µatm) | CO2f (mmol C/m²/d) | DIC (mmol/m³) | Alk (mmol/m³) |
|---|---|---|---|---|---|---|
| 2140-01 | 284.3 | 253.4 | 260.1 | -32.2 | 2006 | 2355 |
| 2140-06 | 284.3 | 198.6 | 205.3 | -25.7 | 1950 | 2355 |
| 2140-12 | 284.3 | 173.2 | 179.8 | -22.4 | 1917 | 2354 |
| 2141-06 | 284.3 | 160.9 | 167.6 | -19.8 | 1896 | 2353 |
| 2141-12 | 284.3 | 151.3 | 157.9 | -19.9 | 1882 | 2352 |

- Atmospheric CO2 arrives correctly (284.25 ppm at every node).
- Surface DIC falls by 124 mmol/m³ in 24 months, alkalinity stays flat: the ocean is losing carbon, not being diluted.
- O2 flux is strongly negative too (-91 mmol O/m²/d in 2140, -26 in 2141).
- pH rises from 8.29 to 8.37.

## The decisive observation

`dpCO2s` is documented as "oceanic pCO2 minus atmospheric pCO2". It equals `pCO2s + 6.7` in every month. So the atmospheric pCO2 used by the flux is **-6.7 µatm**, not about +278.

REcoM (mocsy, `x2pCO2atm`) converts mole fraction to partial pressure as

    pCO2atm = xCO2 * (Patm - pH2O)

With Patm = 0 this is `-xCO2 * pH2O`. Recomputed from the model's own SST and xCO2 (Weiss and Price vapour pressure), January 2141:

- implied by the output (pCO2s - dpCO2s): mean -6.70 µatm, range -12.7 to -1.9
- predicted with Patm = 0: mean -6.59 µatm
- predicted with Patm = 1 atm: mean +277.7 µatm
- node-by-node correlation of the implied value with the Patm = 0 prediction: 0.996

## Cause

`recom_main.F90` picks the sea-level pressure at compile time:

    #if defined(__oasis)
        Loc_slp = pa2atm          ! constant 1 atm, a stopgap written for AWI-ESM2
    #else
        Loc_slp = press_air(n)    ! FESOM's forcing field
    #endif

- FESOM's main target is compiled with `__oasis` and `__oifs`. REcoM is a separate library; its `CMakeLists.txt` passes only `__oifs`, `__recom`, `__usetp`. Checked in `build/src/recom/CMakeFiles/recom.dir/flags.make`.
- So REcoM takes the `#else` branch.
- `press_air` is allocated as zero in `gen_forcing_init.F90` and filled only on the standalone forcing path of `gen_forcing_couple.F90`. In a coupled run nothing writes to it.
- `Loc_slp = 0` goes into `REcoM_Forcing`: `Patm = Loc_slp / Pa2atm = 0` for the CO2 flux, `ppo = 0` for the O2 flux (`o2sat = o2sat_1atm * ppo`, so the saturation is zero), and `Patm_depth = 0` in `recom_sms`.

Winds are not affected. OpenIFS sends 10 m winds and `gen_forcing_couple.F90` puts them into `u_wind`/`v_wind`, which is what REcoM's `#else` wind branch reads.

## Scope

- Every OpenIFS-coupled REcoM build since the coupling went into FESOM main (#806), on any machine.
- In the tagged code: FESOM `2.8.1` pins REcoM f71285c, which has it. `v3.5.0-cc` builds and runs, with wrong carbon and oxygen cycles.
- Physics-only and ice-sheet variants are not affected (REcoM is not compiled).
- Ocean physics of the carbon-cycle variant is not affected in concentration-driven mode: nothing feeds back from REcoM. In an emission-driven run the wrong flux would reach the atmosphere.

## Fix (drafted)

REcoM branch `fix/oifs-coupled-sea-level-pressure`, one commit on REcoM main (ec405f0), local clone at `/albedo/scratch/user/jstreffi/model_codes_v350/recom_fix`; patch in `notes/recom_oifs_sea_level_pressure.patch`.

- `#if defined(__oasis) || defined(__oifs)` for the pressure choice.
- `pa2atm` added to the `use recom_config, only:` list of `recom`. The first build failed without it ("This name does not have a type"): the `__oasis` branch had never been compiled in this configuration.

Why not pass `__oasis` to the REcoM library instead: it would also switch the wind speed to the AWI-ESM2 stopgap that derives 10 m wind from stress, which is worse than the real winds; and it would have hit the same missing import.

## Ideas and open questions

1. **Constant 1 atm is a stopgap.** The real sea-level pressure varies by a few percent (about 0.97 to 1.03 atm in the monthly mean, lower in storm tracks and over the Southern Ocean), and pCO2atm scales with it. The proper fix is for OpenIFS to send mean sea-level pressure and for FESOM to fill `press_air` in `gen_forcing_couple.F90`. That needs a new coupling field in OpenIFS, the namcouple and esm_tools.
2. **Same field, other users.** FESOM's transient tracers (CFC, SF6, 14C, 39Ar in `oce_ale_tracer.F90`) read `press_air` directly. Switched on in a coupled run they would see zero pressure as well. `oce_ale_vel_rhs.F90` uses it for the inverse-barometer term, which is presumably off when coupled; not checked.
3. **A guard would have caught this.** A one-time check in REcoM that `Loc_slp` is within, say, 0.8 to 1.1 atm, aborting otherwise, costs nothing. So would a check that `dpCO2s` is not systematically positive.
4. **How wrong is the existing state.** After 2 years the surface has lost about 6 % of its DIC and the flux is still -20 mmol C/m²/d (about 30 Pg C per year if held globally). Any restart written by an affected run should not be used as a BGC initial state; restart REcoM from climatology after the fix.
5. **What the fixed run should look like.** Cold start from GLODAP/WOA under 284 ppm: pCO2s near 280 µatm in the global mean, dpCO2s of mixed sign and small, CO2f of order ±1 mmol C/m²/d in the mean with the usual pattern (tropical outgassing, high-latitude uptake). The one-month test PICAL_cc_slpfix checks the first of these. Whether the January value of 253 µatm in the broken run was already depressed within the month, or the initial state itself sits below 280, will show there too.
6. **Is anything else compiled away?** Three more preprocessor switches in `recom_main.F90` depend on `__oasis`/`__oifs` (wind, atmospheric CO2 source, one `!defined(__oifs)` block). Worth a read-through of each for the same build-definition mismatch; only the pressure one has been examined.
7. **Was it right on levante before?** Not checked whether earlier AWI-ESM3 carbon-cycle runs on levante (FESOM branch `fesom2.6_recom_awiesm3_co2_coupling`, before #806) built REcoM with `__oasis`. If they did, this is a regression of the merge; if not, those runs are affected as well.

## Release consequence

A REcoM commit, a FESOM submodule bump, a FESOM tag `2.8.2`, and a `v3.5.1-cc` (or `v3.5.1` for all three, to keep one tag set). `v3.5.0` and `v3.5.0-is` are fine as they are.

## Test

PICAL_cc_slpfix: January 2140, same start as PICAL_cc, FESOM = tag `2.8.1` with REcoM f71285c + the fix (binaries in `fesom-bin/recomslp_c9b3582` of the scratch tree `model_codes_v350/awiesm3-develop-cc`). Compare the January means of pCO2s, dpCO2s, CO2f and O2f with the table above.
