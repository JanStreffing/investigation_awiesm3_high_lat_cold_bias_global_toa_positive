# REcoM in the OpenIFS-coupled build: gas exchange at zero sea-level pressure

State of 2026-10-05, late evening. Diagnosis from output and source. Both fixes have passed a one-month test: constant one atmosphere (PICAL_cc_slpfix) and the sea-level pressure of OpenIFS (PICAL_cc_mslp). Tables: `data/clim/recom_slpfix_january2140.txt`, `data/clim/recom_mslp_january2140.txt`. Draft pull request with both REcoM commits: RECOM-Regulated-Ecosystem-Model/REcoM#34. The OpenIFS and FESOM branches and the esm_tools switch are local only.

## Symptom

First biogeochemistry output of the carbon-cycle variant (PICAL_cc_v350val, 2140-2141, REcoM cold start from climatology, concentration-driven at 284.25 ppm). Global means, area-weighted:

| month | pCO2s [uatm] | dpCO2s [uatm] | CO2f [mmol C/m2/d] | surface DIC [mmol/m3] | surface Alk [mmol/m3] |
|---|---|---|---|---|---|
| 2140-01 | 253.4 | 260.1 | -32.2 | 2006 | 2355 |
| 2140-06 | 198.6 | 205.3 | -25.7 | 1950 | 2355 |
| 2140-12 | 173.2 | 179.8 | -22.4 | 1917 | 2354 |
| 2141-06 | 160.9 | 167.6 | -19.8 | 1896 | 2353 |
| 2141-12 | 151.3 | 157.9 | -19.9 | 1882 | 2352 |

- Surface pCO2 falls by 100 uatm in two years and is still falling. It should sit near 280.
- Surface DIC falls by 124 mmol/m3 with alkalinity flat, so carbon is leaving, not being diluted.
- CO2f is defined as the flux into the water. It is negative at every node, every month.
- Oxygen flux is negative everywhere too (-91 mmol/m2/d in the first year, -26 in the second).
- Surface pH rises from 8.29 to 8.37.

## The giveaway

`dpCO2s` is defined as oceanic minus atmospheric pCO2. In every month it equals `pCO2s + 6.7`. So the atmospheric pCO2 the flux saw was **-6.7 uatm**, while the atmospheric mixing ratio in the same output (`xCO2atm`) is 284.25 ppm at every node.

mocsy converts with `pCO2atm = xCO2 * (Patm - pH2O)`. For `Patm = 0` that is `-xCO2 * pH2O`. Recomputed from the model's own SST for January 2141 (Weiss and Price vapour pressure):

| | global mean [uatm] |
|---|---|
| implied by the output, `pCO2s - dpCO2s` | -6.70 (range -12.7 to -1.9) |
| predicted with Patm = 0 atm | -6.59 |
| predicted with Patm = 1 atm | +277.66 |

Node-by-node correlation of the implied value with the Patm = 0 prediction: 0.996. The spatial pattern is the vapour-pressure pattern, i.e. SST.

## Mechanism

`recom_main.F90`:

```fortran
#if defined(__oasis)
            Loc_slp = pa2atm          ! constant one atmosphere, a stopgap written for AWIESM-2.1
#else
            Loc_slp = press_air(n)
#endif
```

- FESOM's own sources are compiled with `-D__oasis -D__oifs`. REcoM is a separate library (`src/recom/CMakeLists.txt`) and gets `-D__oifs -D__recom -D__usetp` only. Checked in the build's `flags.make`.
- So REcoM takes the `#else` branch and reads `press_air`.
- `press_air` is allocated and zeroed in `gen_forcing_init.F90` and filled only in the standalone forcing path of `gen_forcing_couple.F90`. A coupled FESOM never writes to it.
- `Loc_slp = 0` goes to `REcoM_Forcing`: `Patm = Loc_slp / Pa2atm = 0` for the CO2 flux, `ppo = 0` for oxygen (`o2sat = o2sat_1atm * ppo`, so the saturation concentration is zero), and to `REcoM_sms` as `Patm_depth`.

Winds are not affected: the neighbouring switch falls back to `u_wind`/`v_wind`, and those are received from OpenIFS in the coupled path.

## Why nobody saw it

- PICAL_cc wrote no biogeochemistry output at all (the CMIP7 XIOS templates had no REcoM fields until esm_tools fdde4730e). PICAL_cc_v350val is the first run with it.
- In concentration-driven mode nothing feeds back from REcoM to the physics, so every physical score is clean. PICAL_cc and PICAL_cc_v350val are identical in ocean and ice.
- The `__oasis` branch in REcoM has never been compiled in this configuration: it uses `pa2atm`, which `recom_main` does not import. Passing `__oasis` to the library would not have compiled.

## Scope

- Every OpenIFS-coupled REcoM build since the CO2 coupling went into FESOM main (#806), on any machine. That includes PICAL_cc, PICAL_cc_v350val and the TEST_cc* performance runs.
- It is in the tagged release code: FESOM 2.8.1 with REcoM f71285c. `v3.5.0-cc` builds and runs, but its carbon and oxygen cycles are wrong.
- Not affected: the physics-only and ice-sheet variants (REcoM is not compiled), and the ocean physics of the carbon-cycle variant.
- Would matter far more in emission-driven mode, where the ocean flux feeds the atmospheric CO2 tracer.

## Fix, as drafted

One commit in REcoM, `recom_main.F90` only (patch: `notes/recom_oifs_sea_level_pressure.patch`):

- `#if defined(__oasis) || defined(__oifs)` for the constant-pressure branch
- import `pa2atm` from `recom_config`

Local branches, nothing pushed:

- on REcoM `main` (ec405f0): `fix/oifs-coupled-sea-level-pressure`, 941d65b, in `/albedo/scratch/user/jstreffi/model_codes_v350/recom_fix`
- on the commit FESOM 2.8.1 pins (f71285c): c9b3582, in the candidate tree's submodule; FESOM rebuilt from it, binaries in `awiesm3-develop-cc/fesom-bin/recomslp_c9b3582/` (the tagged ones are kept in `fesom-bin/tag_2.8.1/`)

The first draft did not compile (missing import). The second does.

Deliberately not done: defining `__oasis` for the REcoM library. Besides not compiling, it would switch the wind speed to the wind-from-stress stopgap, which is worse than the real 10 m winds REcoM gets now.

## Test

PICAL_cc_slpfix: January 2140, same start as PICAL_cc, FESOM pinned to the fixed build, everything else as PICAL_cc_v350val. Expected in the January mean:

- `dpCO2s = pCO2s - ~278` instead of `pCO2s + 6.7`
- global-mean pCO2s well above the 253 of the broken build
- CO2f and O2f no longer negative everywhere

One month only shows that the pressure is right. Whether the cold-started state then settles near 280 uatm needs years.

## Open points

- **Constant or real pressure.** One atmosphere everywhere ignores the 2-3 % spatial and seasonal variation of sea-level pressure (lower over the Southern Ocean and the subpolar lows), which maps directly onto pCO2atm and O2 saturation. The proper fix is OpenIFS sending mean sea-level pressure and FESOM filling `press_air`; that needs a new coupling field. The constant is what AWI-ESM2 has used.
- **Is the initial state itself in balance?** January 2140 already has pCO2s = 253 uatm as a monthly mean. Part of that is one month of degassing, but the GLODAP/WOA start may also be off for PI. Only visible once the pressure is fixed.
- **Other users of `press_air`.** FESOM's transient tracers (CFC, SF6; `oce_ale_tracer.F90`) read the same field and would see zero in a coupled run. Not switched on here.
- **`Patm_depth` in `REcoM_sms`.** Zero surface pressure also went into the water-column carbonate chemistry. Small next to the hydrostatic pressure at depth, but it was wrong as well.
- **Existing carbon-cycle runs.** Their REcoM states have lost carbon and oxygen and should not be used as restarts for the biogeochemistry. PICAL_cc lost about 120 mmol/m3 of surface DIC in two years; its five-year end state is further gone.
- **Release.** A REcoM commit, a FESOM submodule bump and a FESOM 2.8.2 tag, and by the freeze rule a v3.5.1 for the carbon-cycle variant; v3.5.0 and v3.5.0-is are unaffected.

## Real pressure (2026-10-05, evening)

- OpenIFS `feat/fesom-mslp` 32a073b: sends `A_MSLP` to FESOM under `ECE_CPL_MSLP`.
- FESOM `feat/oifs-mslp-to-recom` 56371fb7: `run_config use_atm_mslp`, receives `mslp_oce` into `press_air`. The momentum-equation pressure term stays under `l_mslp`.
- REcoM 64aaa47 (79f7ae3 on the pinned f71285c): `press_air` where positive, else one atmosphere.
- esm_tools, uncommitted: `general.with_mslp_oce_coupling` (default false) sets both switches and appends `mslp_oce <--gauswgt_i-- A_MSLP` to `rst_co2_ao.nc`.
- Coupler start state with `A_MSLP` = 101325 Pa: `input/oasis_ini/cc_mslp_from_gmhemi1800_2139/`.

PICAL_cc_mslp, January 2140 (`scripts/analysis/recom_mslp_check.py`): the pressure recovered from the output follows the OpenIFS msl (r = 0.997 node by node, rms 0.0009 atm; 0.979 atm over 78-60S, 0.997 global ocean). Against the constant: atmospheric pCO2 -5.9 uatm over 78-60S, -2.6 over 60-40S, -0.8 global. SST and SSS identical.
