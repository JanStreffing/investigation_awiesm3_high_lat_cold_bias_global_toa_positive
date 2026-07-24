# LPJ-GUESS boreal forest-cover bug and the NH high-latitude cold bias — investigation & fixes

*AWI-ESM3 v3.4 (OpenIFS TCO95 + FESOM2 CORE3 + LPJ-GUESS 4.1.2). Companion to the tuning progress report;
all data/figures in `/work/bb1469/a270092/eval/`. Source tree:
`/work/ab0995/a270270/model_codes/awiesm3-develop/lpj_guess/` (identical to `/work/ab0246/a270092/model_codes/awiesm3-v3.4.2/lpj_guess/`).*

## 1. Symptom (measured)
- Coupled runs have a severe **cold-season 2 m cold bias over Siberia / N. Canada** (vs ERA5). Seasonal CMPI `tas`
  maps: bias is **DJF ≈ MAM (deep, saturated < −5 K), collapses to neutral/warm in JJA, moderate in SON** — i.e.
  present while snow is on the ground, gone once it melts (consistent with a snow/vegetation-albedo mechanism),
  but the strong **DJF** part cannot be a shortwave-albedo effect (no boreal sunlight in midwinter) so a second
  cold-season driver (snow insulation / stable BL / longwave) also contributes.
- Boreal high-vegetation cover **cvh ≈ 0.05 with LPJG vs ≈ 0.44 without LPJG** (HTESSEL + prescribed satellite
  veg, AWI-CM3 `TUNE42PI_FES27`), same TCO95 atmosphere, yr1370–79. Per-tile LAI is fine (≈2.5, ≥ satellite).
  **The problem is the forest COVER FRACTION, not leaf density.**
- It is **inherited from the LPJG spin-up**: coupled runs init LPJG from the spin-up's yr-3850 state; baseline
  yr-1350 boreal veg is identical to the spin-up end-state. Not a fast coupled runaway.

## 2. Coupling chain — where cvh comes from (source-verified)
LPJG **computes and sends** the high-veg fraction to OpenIFS (it is *not* prescribed by IFS):
- `namcouple`: `GUE_LLAI:GUE_HLAI:GUE_FRAL:GUE_FRAH:GUE_TYPL:GUE_TYPH → LAILVeg:LAIHVeg:FracLVeg:FracHVeg:TypeLVeg:TypeHVeg` (EXPORTED).
- `GUE_FRAH → FracHVeg` becomes IFS **cvh**; defined `OASIS_Out` in `framework/OasisCoupler.cpp:294`.
- Computed in `framework/framework.cpp::computeDailyLAIandFPCforIFS` (lines 1077–1388). The per-tree-individual
  contribution (framework.cpp:1152–1156):
  ```cpp
  double current_fpc = indiv.crownarea * indiv.densindiv * (1.0 - lambertbeer(indiv.lai_indiv_today()));
  if (iftreefracca && indiv.pft.lifeform == TREE)      // <-- our case
      current_fpc = indiv.crownarea * indiv.densindiv;  // pure crown-area cover, no LAI term
  ```
  summed over tree individuals, area-weighted over patches/stands, capped at `natural_frac` and 1.0
  (framework.cpp:1171,1310,1386). `IFStypehigh` (→ tvh) is a dominance test; over our boreal cells tvh≈3.9
  (evergreen needleleaf) — the *type* is correct, only the *area* is wrong.
- **`iftreefracca = 1`** in the AWI `global.ins` (line 151) → cvh = Σ crownarea·density. (EC-Earth's `ecearth.ins`
  uses `iftreefracca 0`, the FPC+LAI formula.) The much-cited `LOWTOHIGHCUTOFF_FPC = 0.05` constant is **dead
  code** (declared, never used); `natural_frac` = 1 over PI land (no land-use, `file_lunatural ""`). So nothing
  clamps cvh to 0.05 — it is genuinely low **crown area × density**, i.e. boreal trees are far too sparse.

**Conclusion:** cvh≈0.05 ⇒ boreal tree crown cover is genuinely ~5 % (should be ~50 %). This is an
**establishment/density problem**: too few boreal trees establish and survive, so their summed crown area is tiny.

## 3. What gates boreal establishment (the tunable knobs) — `global.ins`
Boreal needleleaf PFTs (verified in `data/ins/global.ins`):
| PFT | tcmin_surv | tcmin_est | tcmax_est | twmin_est | gdd5min_est | notes |
|---|---|---|---|---|---|---|
| BNE (evergreen needleleaf) | −31 | −30 | **−1** | **5** | **500** | shade-intolerant group |
| BINE (interm. shade-tol) | −31 | −30 | **−1** | **5** | **500** | |
| BNS (summergreen needleleaf, larch) | −1000 (none) | −1000 (none) | −2 | −1000 (none) | **350** | greff_min 0.04 |
Global flag: `iftreefracca 1`.

The binding constraints over cold, short-season central Siberia / N. Canada interiors are plausibly
**`gdd5min_est 500`** (growing season too short → BNE/BINE cannot establish) and **`twmin_est 5`** (warmest-month
minimum). Where BNE/BINE are excluded, only BNS (larch) or C3 grass remain; if BNS is also marginal, grass wins →
crown cover collapses → cvh→0.05. This is the classic LPJG "boreal trees replaced by grass under a (cold-biased)
driving climate" failure, self-consistent with the inherited spin-up state.

## 4. Ranked, actionable fixes (respecting the "sunlight & growing-season only" constraint)
1. **Relax BNE/BINE establishment limits (cheapest, most direct).** Lower `gdd5min_est` (500 → ~350, matching
   BNS) and `twmin_est` (5 → ~2–3) for BNE/BINE in `global.ins`. Expectation: boreal trees establish over more of
   Siberia/N. Canada → crown cover (cvh) rises → snow masked → warms spring/shoulder seasons. Cheap to test: a
   short offline LPJG re-spin (the standalone `lpjg-spinup` setup) and read `cvh`/FracHVeg before any coupled run.
   Risk: over-greening / trees where tundra belongs; check against satellite tree-cover.
2. **Re-spin LPJG, don't just restart.** Because the bias is inherited from the yr-3850 state, any establishment
   fix must be applied in a **new LPJG spin-up**; restarting the coupled run from the old state will not help.
3. **`iftreefracca` / crown-area allometry.** With `iftreefracca 1`, cvh = crownarea·density; if crown allometry
   under-sizes boreal crowns, even established trees give low cover. Test `iftreefracca 0` (FPC+LAI, EC-Earth
   setting) as a sensitivity — it changes how sparse canopies map to cover.
4. **Zero-/low-LAI stem snow-masking albedo (OIFS side).** In OpenIFS the snow albedo is a cvh-weighted blend of
   a bright exposed-snow tile and a dark forest-masked tile (`surf/module/surfbc_ctl_mod.F90:389–391` &
   `:430–431`: `PFRTI(5)=ZCVS·(1−cvh)`, `PFRTI(7)=ZCVS·cvh`), with the masked-snow albedo `ZADTI7` a flat
   per-forest-type constant `RALB_SNOW_FOREST` (`surf/module/surfrad_ctl_mod.F90:650–655`; evergreen needleleaf
   = [0.31,0.24]; table set in `susrad_mod.F90:110–125`). Two options: (a) make `ZADTI7` LAI-dependent so leafless
   trees still darken snow — but this needs high-veg LAI threaded into `SURFRAD_CTL` (not currently passed); (b)
   simpler — put a **woody stem-area floor on the masking fraction** in `surfbc_ctl_mod.F90` so tile-7 keeps a
   minimum weight from stems even at low LAI. *Caveat:* because the blend weight is cvh itself, with cvh≈0.05 the
   masked tile gets only 5 % weight — so a stem-albedo tweak mostly helps *after* fix (1) raises cvh, and cannot
   fix the DJF part. (AWI already edits this routine — see `surfrad_ctl_mod.F90:141`, Streffing 2023 EC-Earth
   albedo scaling.)
5. **Don't expect a full fix of the DJF cold.** Spring/shoulder cold is albedo-addressable via 1–4; the co-equal
   DJF cold (no sunlight) points to snow/longwave/boundary-layer and is out of scope for the vegetation lever.

## 5. Prior work
`ecearth.ins` (EC-Earth's LPJG config) is present and differs (`iftreefracca 0`) — worth diffing its boreal
settings. No explicit "Laszlo" experiment directory was located in the searched trees; recommend confirming with
him directly whether establishment-limit tests were already run.

## 6. Bottom line
LPJG under-establishes boreal forest → sends cvh≈0.05 (vs ~0.44 observed) to OpenIFS → snow unmasked →
cold-season (esp. spring) NH high-latitude cold bias, shared by all tuning runs and inherited from the spin-up.
The highest-value, low-cost fix is to **relax BNE/BINE establishment limits (`gdd5min_est`, `twmin_est`) and re-do
the LPJG spin-up**, verifying that boreal `cvh` recovers toward the satellite/HTESSEL ~0.44. This is the
physically-grounded *warming lever* that can offset the imbalance-tuning cold bias.

## 7. 2026-07-24 coupled follow-up — what changed after a 50-year run

The newer coupled follow-up from the CRUNCEP3-initialized branch (alias "080a")
changes the interpretation in one important way: the boreal problem is **not adequately described as a purely climate-free,
competition-only issue once the coupled system is allowed to evolve for 50 years**.

Key new findings from the finished 50-year run:

- NH45+ `TREEFPC` still declines strongly (`0.336 -> 0.215`).
- The collapse is strongest over **Siberia / East Siberia**, where `TREEFPC` and
  `AGDD5` both decrease strongly by year 50.
- At the same time, NH45+ aggregate `AGDD5` increases, and some other regions
  (e.g. Scandinavia / Canada) do not show the same thermal decline.
- The final-decade seasonal temperature bias map against CRUNCEP3 shows a
  **warm-ocean / cold-land split**, not a simple hemispheric cold drift.

Updated interpretation:

- The older decisive test still supports the claim that the **spin-up state and
  LPJG competition structure matter a lot**.
- But the 50-year coupled follow-up shows that **regional coupled climate also
  matters**, especially over boreal continental land.
- Therefore the current best statement is:
  - the boreal cold bias has a **structural LPJG component**,
  - but the coupled system can still generate a **regional land-cooling pattern**
    that amplifies the Siberian forest loss.

Practical consequence:

- Do **not** treat vegetation-only retuning as sufficient.
- The next tuning path should combine:
  1. a better coupled physical branch (`06T` / `06V`) to improve radiation
     balance without worsening East Siberian land climate, and only then
  2. LPJ competition retuning (`07A`, then `07A+07C`) to retain boreal trees.

See also:

- `notes/COUPLED_CRUNCEP3_INIT_BRANCH_FOLLOWUP_2026-07-24.md`

*Verified-from-source: coupling chain (§2), iftreefracca=1, establishment limits (§3). Inferred: that
gdd5min_est/twmin_est are the binding constraints (consistent with the maps but not yet confirmed by an
`est_limits.out` analysis — a good next diagnostic).*
