# AWI‑ESM3 (v3.4, CORE3 / TCO95) PI‑control tuning campaign — `Tuning_test_06*`

## Addendum (2026-07-28): new-sea-ice branch `Tuning_test_09*`

This report now also tracks the new-sea-ice continuation branch:

- `09A`: baseline coupled run with the new sea-ice scheme
- `09B`: new sea-ice + tuned coupled run (`06T` style stack)
- `09C`: new sea-ice coupled `06V` variant

### Run status and window used

All three runs reached the target 30-year segment (`1350-1379`) and are available for
consistent last-10-year diagnostics (`1370-1379`):

- `09A`: `Tuning_test_09A_lpjguess_Baseline_coupled_fromCRUNCEP_newSeaIce`
- `09B`: `Tuning_test_09B_06T_1hCPL_MOSPP_KPPLOW_CRUNCEPinit_newSeaIce`
- `09C`: `Tuning_test_09C_06V_CRUNCEPinit_newSeaIce`

### Updated comparison figures

Generated with common plotting scripts in:
`/work/ab0995/a270270/analysis/LR_offline_LPJ_GUESS_tunning_new/scripts/`

- T2m (ocean + land):
	`/work/ab0995/a270270/analysis/LR_offline_LPJ_GUESS_tunning_new/plots/diff_from_amip_and_era5/t2m_bias_vs_cruncep3_spinup_09A_09B_09C_nh_last10yrmean.png`
- TREEFPC:
	`/work/ab0995/a270270/analysis/LR_offline_LPJ_GUESS_tunning_new/plots/diff_from_amip_and_era5/treefpc_bias_vs_cruncep3_spinup_09A_09B_09C_nh_last10yrmean.png`

Both figures use:

- spinup reference for vegetation (`CRUNCEPcalibrated_v3`, last 10 years in spinup file)
- last-10-year means for coupled runs from model years `1370-1379`

### `09B` tuning block documented on figures (original -> tuned)

- `oasis3mct.time_step`: `7200 s -> 3600 s`
- `use_momix`: `.false. -> .true.`
- `kpp_av0`: `0.005 -> 0.003`
- `kpp_kv0`: `0.005 -> 0.003`
- `kpp_avbckg`: `1e-4 -> 5e-5`
- `kpp_kvbckg`: `1e-5 -> 5e-6`
- `pndaspect`: `0.8 -> 1.3`
- `rfracmax`: `1.0 -> 0.75`
- `albpnd`: `0.2 -> 0.28`

### Interpretation update

The active branch comparison is now `09A` vs `09B` vs `09C` under the same new-sea-ice
framework and year window. This supersedes using only `09A/09B` for current
decision-making and should be the default set in follow-up discussions.

**Progress report** · prepared 2026‑06‑30
**Experiments by:** a270270 (`/work/bb1469/a270270/runtime/awiesm3-v3.4/`)
**Evaluation:** release_evaluation_tool2 (`/work/ab0246/a270092/software/release_evaluation_tool2`)
**All figures & data:** `output/Tuning_test_06_overview/`

---

## 1. Goal

Reduce the **net top‑of‑atmosphere (TOA) radiative imbalance** of the AWI‑ESM3 pre‑industrial control so the model reaches equilibrium faster, **without**

* degrading the **CMPI** score (CMIP6‑normalised multi‑variable performance index, part4), and
* inflating the **2 m‑temperature RMSE vs ERA5** (part8 bias map).

The baseline pre‑industrial control still drifts with a **positive TOA imbalance of ≈ +1.5 W/m²** (energy flowing into the system), which is the quantity being tuned down.

## 2. What was run

22 experiments, each a 30‑year segment (model years **1350–1379**) of the CORE3 / TCO95 configuration (FESOM2 CORE3 mesh, 211 567 nodes; OpenIFS TCO95 / A096; LPJ‑GUESS). All share one **Baseline** and vary one or a few parameters.

> **Important structural note.** Nearly every non‑baseline run also carries a *common melt‑pond shift* relative to the true Baseline: `albpnd 0.2→0.28`, `rfracmax 1.0→0.75`, `pndaspect 0.8→1.3` (FESOM `namelist.ice &meltpond`). The A/B/C runs build this up incrementally; from D onward it is the de‑facto starting point on top of which each run adds its own change. So differences quoted below are *each run's distinguishing change*, on top of that common pond shift (except where noted).

### Parameter catalogue

| Run | Distinguishing change (file / namelist) | Baseline → New | Process |
|---|---|---|---|
| **06_Baseline** | — | — | reference |
| 06A | `albpnd` (`namelist.ice &meltpond`) | 0.2 → 0.28 | sea‑ice melt‑pond albedo |
| 06B | + `rfracmax` | 1.0 → 0.75 | pond meltwater retention |
| 06C | + `pndaspect` | 0.8 → 1.3 | pond depth:area (= "common pond" combo) |
| 06D *HRlike* | (only the common pond shift) | — | ⚠ no other namelist diff found (see §7) |
| 06G | `albpnd / rfracmax / pndaspect` | 0.35 / 0.6 / 1.6 | "stronger meltpond" |
| 06L | `albpnd / rfracmax / pndaspect` | 0.4 / 0.55 / 1.8 | "extra strong meltpond" |
| 06H | strong pond (like G) + `RVICE` (`&NAMCLDP`) | 0.16 → 0.18 | meltpond + cloud‑ice fall speed |
| 06Q | `albw`/`albocn` | 0.1 → 0.045 | open‑water / ocean albedo |
| 06R | `albw`/`albocn` | 0.1 → 0.075 | open‑water / ocean albedo |
| 06S | `whichevp` (`&ice_dyn`) | 1 → 0 | sea‑ice rheology (mEVP→EVP) |
| 06E | `RVICE` + `GGAUSSB` (`&NAMGWWMS`) | 0.18 ; −0.5→−0.6 | cloud‑ice + non‑oro gravity‑wave |
| 06F | `LRDALB` (`NALBEDOSCHEME=3`) | F → T | diagnosed (MODIS‑type) surface albedo |
| 06P | `ENTSTPC3` (`&NAMCUMF`) | (absent) → 1 | deep‑convection entrainment |
| 06N | `use_momix` (`namelist.tra`) | F → T | near‑surface ocean mixing ("mospp") |
| 06M | OASIS coupling period + LAG | 7200 s → 3600 s | atmosphere↔ocean coupling 2h→1h ("1hcpl") |
| 06O | 1hcpl + mospp | — | coupling + ocean mixing |
| 06O4 | O + `momix_kv` | 0.01 → 0.012 | + larger Monin‑Obukhov mixing |
| 06O5 | O + `momix_kv` | 0.01 → 0.02 | + even larger mixing |
| 06T | O + KPP‑low (`kpp_av0/kv0` 0.005→0.003, bckg ↓10×) | — | + weaker KPP vertical mixing |
| 06U | T + open‑water albedo 0.075 | — | combo |
| 06V | T + `ENTSTPC3=1` | — | combo (most aggressive) |

## 3. Headline result — the imbalance ↔ bias trade‑off

> **Figure:** `plots/summary_scatter_official.png` (x = part2 TOA imbalance, y = part8 ERA5 RMSD; dotted lines = baseline; lower‑left is better).
> **Figure:** `plots/summary_TOA_timeseries.png` (annual TOA imbalance, all 22 runs).

Every change that meaningfully lowers the TOA imbalance does so by **cooling the model**, which simultaneously **worsens the cold bias against ERA5**. The two goals pull against each other along a clear front:

* **Sea‑ice / albedo / meltpond family (A,B,C,D,G,L,Q,R,S, + E,F,H,P):** modest imbalance reduction (5–20 %), small RMSD penalty (+0.1 to +0.4 K). These move *gently* down‑left.
* **Coupling + ocean‑mixing family (M,O,O4,O5,T,U,V):** large imbalance reduction (30–57 %) but a **large cold bias** penalty (RMSD +0.6 to +0.95 K, mean bias down to −2.1 K).

### Master metric table (sorted by TOA imbalance, best first)

| Run | TOA imbalance [W/m²] | vs base | 2m RMSD vs ERA5 [K] | 2m bias [K] | CMPI |
|---|---:|---:|---:|---:|---:|
| 06V_…_entstpc3_1 | **0.655** | −57 % | 2.671 | −2.124 | _pending_ |
| 06T_…_kpplow | 0.952 | −38 % | 2.385 | −1.791 | _pending_ |
| 06O5_…_kv002 | 0.968 | −37 % | 2.511 | −1.904 | — |
| 06O_1hcpl_mospp | 0.978 | −36 % | 2.389 | −1.774 | _pending_ |
| 06O4_…_kv0012 | 1.007 | −34 % | 2.387 | −1.807 | — |
| 06U_…_openwater | 1.020 | −33 % | 2.422 | −1.845 | — |
| 06M_1hcpl | 1.062 | −31 % | 2.315 | −1.658 | — |
| 06P_entstpc3_1 | 1.091 | −29 % | 2.076 | −1.358 | — |
| 06H_combo_g_rvice018 | 1.220 | −20 % | 2.082 | −1.394 | _pending_ |
| 06G_stronger_meltpond | 1.274 | −17 % | 1.990 | −1.274 | — |
| 06L_extra_strong_meltpond | 1.320 | −14 % | 2.050 | −1.332 | — |
| 06E_rvice018_ggaussb | 1.335 | −13 % | 1.969 | −1.231 | — |
| 06D_HRlike | 1.377 | −10 % | 1.892 | −1.146 | _pending_ |
| 06F_lrdalb_true | 1.378 | −10 % | 1.962 | −1.219 | — |
| 06N_mospp | 1.382 | −10 % | 2.011 | −1.310 | — |
| 06C_albpnd028_rfrac075_pndaspect13 | 1.391 | −9 % | 1.884 | −1.109 | — |
| 06Q_ralbsead0045 | 1.397 | −9 % | 1.891 | −1.154 | — |
| 06S_evp0 | 1.419 | −7 % | 1.935 | −1.191 | — |
| 06R_openwater_albedo0075 | 1.443 | −6 % | 1.860 | −1.127 | — |
| 06B_albpnd028_rfrac075 | 1.447 | −5 % | 1.847 | −1.062 | — |
| 06A_albpnd028 | 1.457 | −5 % | 1.811 | −1.004 | 0.748 |
| **06_Baseline** | **1.531** | — | **1.717** | **−0.939** | **0.745** |

*TOA imbalance = 30‑yr global‑mean net TOA (tsr+ttr) from the tool's `part2_rad_balance.py`. RMSD/bias = `part8_t2m_vs_era5.py` methodology (cos‑lat‑weighted, model − ERA5, model yrs 1354–1379); cross‑checked against the tool's on‑figure values (baseline 1.717 vs 1.721; 06V 2.671 vs 2.678). CMPI from `part4_cmpi.py` (< 1 = better than CMIP6 ensemble‑mean error; lower better).*

## 4. Radiative imbalance (part2)

Per‑run `radiation_budget.png` plots are in `output/<run>/`. Key points:

* Baseline net TOA imbalance ≈ **+1.53 W/m²** (30‑yr mean), still clearly positive → ongoing warm drift.
* The **surface** net flux tracks TOA almost exactly (TOA − SFC ≈ −0.09 W/m² in all runs), i.e. the imbalance is genuinely a TOA radiation problem, not a spurious surface‑flux leak.
* Strongest reductions come from **1‑hour coupling + near‑surface ocean mixing (mospp) + weaker KPP**, culminating in 06V at **+0.66 W/m²** (≈ 57 % smaller than baseline).
* All runs still trend **upward** over the 30 years (≈ +0.3 to +0.5 W/m²/decade) — none has fully flattened, so 30 years is too short to declare equilibrium; the ranking by mean level is the usable signal.

## 5. 2 m temperature vs ERA5 (part8)

Per‑run `t2m_vs_ERA5.png` maps in `output/<run>/`; candidate comparison in `plots/t2m_bias_multipanel.png`.

* **Baseline** (`output/Tuning_test_06_Baseline/t2m_vs_ERA5.png`): RMSD 1.72 K, bias −0.94 K — cold over NH continents (Siberia, N America, Himalaya, Sahara), warm over the Southern Ocean and tropical upwelling.
* **06V** (`output/…_entstpc3_1/t2m_vs_ERA5.png`): RMSD 2.68 K, bias −2.12 K — the map is **cold almost everywhere**. The imbalance was removed by globally cooling the planet, which over‑corrects the Southern‑Ocean warm bias into a pervasive cold bias.
* **Caveat:** this is a **PI‑control vs present‑day ERA5** comparison, so some global cold offset (≈ 0.5–0.7 K) is *physically expected*. The baseline's −0.94 K is already partly this; the coupling runs' −1.7 to −2.1 K are clearly **too cold**. The meltpond family (−1.0 to −1.4 K) sits in a more defensible range.

## 6. CMPI (part4) — the tie‑breaker

CMPI integrates 13 variables (siconc, tas, clt, pr, rlut, uas, vas, ua@300, zg@500, zos, mlotst, thetao, so) over regions/seasons, normalised against 30 CMIP6 models; **< 1 means better than the CMIP6 ensemble‑mean error.**

Run on a representative subset (Baseline, 06A, 06D, 06H, 06O, 06T, 06V) spanning the trade‑off front.

| Run | CMPI |
|---|---:|
| 06_Baseline | 0.745 |
| 06A_albpnd028 | 0.748 |
| 06D_HRlike | _pending_ |
| 06H_combo_g_rvice018 | _pending_ |
| 06O_1hcpl_mospp | _pending_ |
| 06T_1hcpl_mospp_kpplow | _pending_ |
| 06V_…_entstpc3_1 | _pending_ |

*(CMPI job 25962939 still running; this section is updated as each run completes. Heatmaps land at `output/<run>/cmpi.png`.)*

The baseline is already a good model (CMPI 0.745). The key question CMPI answers: do the strong‑cooling runs (O/T/V) keep CMPI acceptable, or does the cold bias drag it up? Their ERA5 T2m RMSD roughly doubled the bias, so the expectation is a CMPI penalty — but CMPI's multi‑variable normalisation may absorb part of it.

## 7. Caveats / things to flag to the experimenter

1. **06D "HRlike" carries no HR‑specific namelist change** — the only diff vs Baseline is the common pond shift (identical to 06C). If "HRlike" was meant to change resolution/timestep/viscosity, that is **not present** in any namelist diffed; its near‑identical metrics to 06C confirm this.
2. **Filename convention.** These runs write monthly remapped fields as `atm_remapped_1m_<var>_YYYY-YYYY.nc` (no `_1m_` frequency infix). The tool's `part2`, `part8` and the CMPI preprocessing only knew the `_1m_`‑infix name, so I added **backward‑compatible fallbacks** to `scripts/part2_rad_balance.py`, `scripts/part8_t2m_vs_era5.py`, and `preprocessing_examples/preprocess_AWI-CM3-XIOS.sh`. Existing configs are unaffected.
3. **Accumulation period.** These XIOS runs accumulate fluxes over **3600 s** (units J m⁻², `interval_operation = 3600 s`), not the 21600 s of the old HR config. Configs use `accumulation_period = 3600` (verified: tsr → ≈ 240 W/m², physical). Using 21600 would have under‑stated all fluxes by 6×.
4. Imbalance is read over the full 30‑yr segment; none of the runs is at equilibrium yet, so absolute numbers will keep drifting — treat them as a **ranking**, not final equilibrium values.

## 8. Recommendation

* The **coupling/mixing runs (O/T/V) are the only way found so far to cut the imbalance substantially (>30 %)**, but they overshoot into a strong global cold bias. Whether they are acceptable hinges on their **CMPI** (pending) and on how much PI‑vs‑present cold offset you are willing to allow.
* If CMPI of O/T/V degrades, the **best‑balanced candidates are 06P (entstpc3, −29 % imbalance, RMSD +0.36) and 06H (meltpond+rvice, −20 %, +0.37)** — meaningful imbalance reduction at a fraction of the cold‑bias cost.
* Consider **combining a moderate coupling/mixing setting with a *warming* lever** (e.g. lower cloud/SW reflection) to recover the cooling‑induced bias while keeping the imbalance gain — the current combos only stack *cooling* levers.

---

### Reproduce
```
source ~/loadconda.sh && conda activate reval
cd /work/ab0246/a270092/software/release_evaluation_tool2
python scripts/part2_rad_balance.py output/Tuning_test_06_overview/configs/<run>.py
python scripts/part8_t2m_vs_era5.py  output/Tuning_test_06_overview/configs/<run>.py
python scripts/part4_cmpi.py         output/Tuning_test_06_overview/configs/<run>.py
```
Per‑run configs: `output/Tuning_test_06_overview/configs/`. Master data: `output/Tuning_test_06_overview/data/MASTER_metrics.csv`.
