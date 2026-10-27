# CO2 emissions for the emission-driven historical run: plan and status

Written 2026-10-08. Purpose: AWI-ESM3 needs anthropogenic CO2 emissions for `esm-hist` (and later
`esm-scen7-*`). This note records the decision on how they enter the model and the work it takes.
A first implementation exists since 2026-10-08, see "Status" at the end; nothing is pushed.

## Decision

- The emissions go into the **OpenIFS CO2 tracer**, as in EC-Earth4.
- They arrive **through the coupler** from a small program, the **forcing reader**, not through a GRIB file read by
  OpenIFS. EC-Earth4 uses the file route (`ICMCL<exp>INIT_COMPO` with `LMCC_COMPO`); we do not.
- The forcing reader is a **Fortran derivative of the AMIP forcing reader**. "AMIP" implies no
  ocean, which is wrong for a coupled run, and "data coupler" is too wide: it reads and sends in one
  direction. The name says what it does and leaves room for other prescribed fields later
  (peat-fire CO2, deposition, other emissions, SST and sea ice for atmosphere-only runs).
  Names: repository `forcing_reader`, executable and esm_tools component `freader`, OASIS grid
  `FR05`, sent fields `FR_CO2FF` and `FR_CO2AIR01` to `25`.
- **Aircraft emissions keep their 25 height levels**: 25 fields plus 1 surface field, exchanged once
  a day. The number of fields is set by the data, not by the model: with L91, 51 model levels
  (41 to 91) lie in the file's range of 0 to 15.25 km, and every 0.61 km band contains at least one
  (one or two above 3 km, nine in the lowest band). OpenIFS spreads each band's field over the model
  levels inside it. Merging to about 5 bands would be hard to tell apart in any result for a
  long-lived gas, but 26 fields once a day cost little, so all 25 are kept.

Why the coupler and not the file:

- OpenIFS already receives CO2 fluxes this way in our `-cc` setup (one field from REcoM, three from
  LPJ-GUESS). A further field uses a mechanism that runs today.
- The input stays the NetCDF as published. OASIS does the conservative remapping; there is no GRIB
  converter and no regridding code to write.
- No file has to be cut for each leg, and a new resolution or scenario needs no preprocessing.
- With `LMCC_COMPO` on, OpenIFS expects every time-varying emission field in the file
  (`sumcclag.F90`), including our two coupled sectors. That problem does not arise.

What it costs: an OpenIFS code change that has to go to the consortium branch and into a new tag,
one more executable in the coupled launch, and a route that differs from EC-Earth4's. EC-Earth4
prototyped the coupler route in 2024 (ec-earth-data-coupler issue 2, field `A_CO2_EMIS`) and later
moved to the file; the reason is not recorded on GitLab. **Ask Etienne Tourigny before building.**

## Input data (downloaded, in the pool)

`/albedo/pool/oifs-48r1/cmip7-data/emissions/`

| folder | content | status |
|---|---|---|
| `CEDS-CMIP-2025-04-18/` | historical, 1750-2023: `CO2_em_anthro` (6 files) and `CO2_em_AIR_anthro` (6 files), plus cell area | complete, SHA-256 verified against ESGF |
| `IIASA-IAMC-h/` | scenario `h`, 2022-2100: surface (`-1-1-1`) and aircraft (`-1-1-2`) | complete, SHA-256 verified |
| `IIASA-IAMC-vl/` | scenario `vl`, the same two files | surface file complete and verified; aircraft file still downloading on 2026-10-08 (434 of 1183 MB) |

Facts read from the files:

- Grid 0.5 degree (720 x 360), monthly, units kg m-2 s-1 of CO2, calendar `365_day`.
- Surface file: 8 sectors (agriculture, energy, industrial, transportation, residential and
  commercial, solvents, waste, international shipping). Fossil and industrial CO2 only.
- Aircraft file: 25 levels of 0.61 km, centres 0.305 to 14.945 km. Zero before 1900, about 2.5 % of
  the total in 2023.
- Scenario surface files have 10 sectors (two for carbon removal, which can be negative) and hold
  values for selected years only; months have to be interpolated between those years.

Not used: biomass-burning CO2 (`DRES-CMIP-BB4CMIP7`). LPJ-GUESS simulates fire and land use, so
adding it would count deforestation carbon twice. Peat fires are the open exception, see below.

## What OpenIFS can already do

Read in the `48r1v6` source:

- A 2D emission entry (`NAMCOMPO_EMIS`) has a `VERTICAL_PROFILE_TYPE`: `Surface`, `HeightRange`,
  `AltitudeRange`, `LevelRange`, `HeightMap`, `AltitudeMap`, `GFAS`
  (`compo_injection_profile.F90`). So a 2D field can be injected into an altitude band by namelist.
- Up to 500 2D emission fields are allowed (`NPEMIS2D`), so 26 are no problem.
- True 3D emission fields exist too (`APPLY_3D_EMISSION`), but nothing fills them from the coupler
  and OASIS does not remap vertically. Not needed for this plan.
- The coupled CO2 fields are received in `ece_recom_get_fluxes.F90` and `ece_lpjg_get_fluxes.F90`
  and written into the tracer's sectors by sector name. The new receive routine mirrors these.

## Work packages

### 1. The reader (Fortran)

- Fork `ec-earth/amip-forcing-ece4` (DKRZ GitLab) under a new name. About 600 lines today, with two
  hard-wired fields (SST, sea ice).
- Read a list of fields from the namelist: file, variable, coupling name.
- Keep the existing options: fixed year (for `esm-piControl`), time interpolation of monthly values.
- Add: several input files in sequence (historical, then scenario), the 365-day calendar of the
  input against the model's leap years, interpolation between the given years of a scenario file.
- Grid 720 x 360 with cell areas, written to the OASIS grid files.

### 2. Input preparation (one script, can live in ocp-tool)

- Sum the 8 (or 10) surface sectors into one field.
- Leave the 25 aircraft levels as they are.
- Write one NetCDF per data source in the form the reader expects.
- Check the annual global totals against the `global_total_emission_*` attributes of the originals.

### 3. OpenIFS

- Declare 26 new coupled input fields and a routine that writes them into emission sectors.
- Namelist: one `NAMCOMPO_EMIS` entry for the surface field in sector `co2apf` (paramId 210069), and
  25 entries for aircraft, each with `VERTICAL_PROFILE_TYPE='AltitudeRange'` and its 0.61 km band.
- Check that LPJ-GUESS's anthropogenic land flux stays in its own sector and is not overwritten.
- Merge request to the consortium's `ecearth` branch, then `awiesm3-develop` and a new tag.

### 4. esm_tools

- A new component for the reader, added to the `-cc` couplings behind a switch for emission-driven
  runs.
- 26 coupling entries, conservative remapping from the 0.5 degree grid to the OpenIFS grid with
  OASIS's global conservation, once a day, with a coupler restart entry for the lag.
- The namelist entries of package 3, generated, not written by hand.

### 5. Tests, in this order

1. Reader alone with OASIS, one month: fields arrive, global integral before and after remapping
   agree.
2. Coupled, one year, emissions fixed at 1850: atmospheric CO2 mass grows by the emitted amount
   (191 Tg CO2 in that year) on top of the land and ocean fluxes.
3. Coupled, a year after 1950: the aircraft part appears at the right altitudes.
4. A run across 2021/2022 with scenario `h`: the reader changes file without a jump.

## Open points

- **Why EC-Earth4 left the coupler route.** One question to Etienne; a blocker there changes this
  plan.
- **Home of the repository.** It is a local git repository on albedo so far
  (`/albedo/work/user/jstreffi/software/forcing_reader`); it needs a place Jan controls before
  esm_master can fetch it.
- **Tracer mass conservation in OpenIFS** is still open upstream (ecearth4 MR 300, oifs48r1.1 MR 24).
  It limits what any emission-driven run is worth, whichever route delivers the emissions.
- **Peat-fire CO2.** Not simulated by LPJ-GUESS. Historical data exist, scenario data do not. The
  EC-Earth thread (ecearth4 issue 261) leans to a constant climatological cycle; no decision.
- **Protocol.** No CMIP7 or C4MIP text was found that prescribes gridded versus global fossil CO2,
  aircraft at altitude, or cement carbonation. Searched through summaries, so "not found", not
  "does not exist".
- **The reader in a five-model launch.** It has only run in the atmosphere-only setup. Process
  placement and the extra remapping weights in the 36-node layout are untested.
- **LPJ-GUESS scenario inputs** are not wired in esm_tools yet; needed for any scenario run,
  independent of this plan.

## Sources

- EC-Earth GitLab: ecearth4 issues 131, 261, 303, 15 and MR 300; ecearth3 issues 1013, 1016;
  ec-earth-data-coupler issue 2.
- ESGF queries of 2026-10-08 for `CEDS-CMIP-2025-04-18` and `IIASA-IAMC-{h,vl}-1-1-1/-1-1-2`.
- input4MIPs dataset overviews: anthropogenic SLCF and CO2 emissions; open biomass burning emissions.
- OpenIFS source at tag `48r1v6`: `updclie_compo.F90`, `sumcclag.F90`, `compo_apply_emissions.F90`,
  `compo_injection_profile.F90`, `yomcompo.F90`, `ecearth/couple_*`.

## Status (2026-10-08)

Built, all local, nothing pushed:

- **Forcing reader**: `/albedo/work/user/jstreffi/software/forcing_reader`, one commit. Reads the
  files as published, so work package 2 (input preparation) is not needed: the sector sum, the
  monthly hold, the interpolation between scenario years and the hand-over from historical to
  scenario files happen in the reader. Its log gives the global total of every month read.
- **OpenIFS**: branch `pr13/datacpl-co2-emissions` off `ecearth` (the branch name predates the
  rename), and the same commit on `awiesm3-develop` in the scratch build tree
  (`model_codes_v350/awiesm3-develop-cc`, branch `test/datacpl`). Switch `ECE_CPL_FREADER_CO2`,
  level count `ECE_FREADER_CO2_NLEV`, routine `ECE_FREADER_GET_FLUXES`, received fields `FCO2FF` and
  `FCO2AIR` (25 levels), tracer sectors `co2ff` and `co2air01` to `co2air25`.
- **esm_tools** (uncommitted): component `configs/components/freader`, switch
  `general.with_co2_emissions`, 26 coupled fields with first-order conservative remapping once a
  day, the 26 `NAMCOMPO_EMIS` entries (`namcompo_emis_co2_emi.nml`).

Found on the way:

- **The atmosphere has no unmasked grid in the coupler files.** `A096` masks land, `R096` masks
  ocean; emissions are on both. An unmasked copy `F096` is needed in `grids.nc`, `masks.nc` and
  `areas.nc`. For the tests it is in a private copy
  (`/albedo/scratch/user/jstreffi/icb_eval/oasis_datacpl_test`); for production ocp-tool should
  write it and it has to go to the pool.
- The reader does not add a node to the 36-node layout: its single task takes a free core.

Not yet done: the twin test (January 2000, with and without emissions), a conservation check of
the remapping, the fixed-year and scenario modes, the merge request, and everything under "Open
points".


## Status (2026-10-27): done for v3.5.1-cc

Renamed to **emissions reader** everywhere (component `emissions_reader`, OASIS `EMISREAD`/`EM05`/`EM_CO2*`,
OpenIFS switch `ECE_CPL_EMIS_CO2`).

- Repository: git.smhi.se/jan.streffing/emissions_reader (private), tag `v1.0`.
- OpenIFS: tag `48r1v7` moved to `76fdce0` (e65d859 + receiving code); also on `awiesm3-develop` (`40fa031`).
- esm_tools `26ea10154`: component in the develop, v3.5.1 and oifsamip cy48 couplings; `general.with_co2_emissions`
  in the runscript; emission-driven runs give LPJ-GUESS the tracer CO2; tracer cold start from `oifs.co2_init_mmr`.
- Pool: `F096` in `oasis/cy48r1/TCO95-CORE3` and `TCO95-AMIP` (old files in `replaced_20261027/`), `AMIP` grid
  in `TCO95-AMIP`, `oifs-48r1/TCO95L91/ICMGGab45INIT_AMIP`, `oifs-48r1/amip-forcing/`.
- Tests: v3.5.1-cc branched from the v3.5.1 spin-up at 2365, one month with emissions (year 2000): stable, tracer
  284.0 -> 284.4 ppm. oifsamip pair with/without emissions: found/emitted 0.996 after one month, the shortfall
  growing with the excess, so inside OpenIFS transport/fixer, not ocean uptake.
- Open: LPJ-GUESS fire spikes up to 15,000 ppm at single African savanna points in January; the tracer's
  transport/fixer leak; aviation 4 % (about 0.1 % of the total); other resolutions need `F<nnn>` in their grid files.
