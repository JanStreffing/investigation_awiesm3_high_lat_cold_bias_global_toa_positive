# AWI-ESM3 v3.5.0: component tags

State of 2026-10-05. All nine tags exist (FESOM `2.8.1` cut by Jan; the other eight pushed the same day as annotated tags), each verified on its remote as a tag only, on the commit listed. Remaining work is the esm_tools side and the validation.

## Rules

- esm_tools must point at tags only, never at branches. The esm_tools key is called `branch:`, but for a release its value is a tag name.
- A tag name must not also exist as a branch in the same repository. `git clone -b` prefers the branch, so the pin would silently move.
- FESOM tags are cut from `main`.
- v3.5.0 is frozen once the esm_tools coupling is committed. A later fix is v3.5.1 with its own coupling file. (The v3.4.2 coupling was bumped four times after it was first pinned.)

## Tags (all created 2026-10-05)

One set for all three variants (physics-only, carbon-cycle, ice-sheet).

| Component | Tag | Remote | Branch | Commit today |
|---|---|---|---|---|
| FESOM | `2.8.1` **(tagged 2026-10-05)** | `https://github.com/FESOM/fesom2.git` | `main` | 44e4f9c0 |
| OpenIFS | `48r1v5` | `https://git.smhi.se/jan.streffing/oifs48r1.1.git` | `awiesm3-develop` | f6d1860d |
| LPJ-GUESS | `4.1.13` | `https://git.smhi.se/jan.streffing/lpjg-4.1.git` | `lpj_guess_awiesm3` | 179fb0d6 |
| OASIS3-MCT | `5.3` | `https://git.smhi.se/jan.streffing/oasis3-mct-5` | `local_combined_fixes` | 898c8b2e |
| XIOS | `2.5.4` | `https://gitlab.dkrz.de/ec-earth/xios-2.5-ece.git` | `main` | 9c33506d |
| Runoff mapper | `1.4` | `https://gitlab.dkrz.de/ec-earth/runoff-mapper.git` | `ismm-freshwater` | 36c888f7 |

Ice-sheet variant only. These are forks that will not merge back, so the tag is the upstream version it is based on plus `-awi.N`:

| Component | Tag | Remote | Branch | Commit today |
|---|---|---|---|---|
| PISM | `1.2.0-awi.1` | `https://github.com/JanStreffing/pism.git` (fork of xihao001/pism, itself a fork of pism/pism) | `crodehacke/detect_holes-awiesm3-is` | 1edf404c |
| dEBM | `1.0-awi.1` | `https://github.com/JanStreffing/dEBM.git` (fork of ukrebska/dEBM) | `esm_tools` | 1889719e |
| ismm | `4.1.8-awi.1` | `https://git.smhi.se/jan.streffing/ecearth4.git` | `cissembel-levante-intel` | 94511fb3 |

Distance from the last tag: FESOM `main` is 140 commits past `2.8.0`, LPJ-GUESS 22 past `4.1.12`, OASIS 6 past `5.2`, XIOS 1 past `2.5.3` (a Levante link flag). The OpenIFS repository `oifs48r1.1` has no tags yet; `48r1v5` exists there neither as tag nor as branch.

Runoff mapper: `ismm-freshwater` is the standard branch (`enthalpy_on_ocean_side`, tagged `1.3`) plus one commit, the EC-Earth4 runoff mapper for the ice-sheet freshwater return. With `NumCoupledRegions = 0` the mapping adds zero fields and multiplies by scaling factors of 1.0, so it should give the same runoff and calving as `1.3`. That is from reading the code: no run has used it yet, it needs the `1.4-ismm` namelist, and its OASIS component name changes from `rnfma` to `RNFMAP`.

PISM: the branch is upstream `v1.2` (= `v1.2.0`) plus 40 commits; it is 19 commits behind `v1.2.1`, so the esm_tools name `pism-github1.2.1` implies the wrong base and should be renamed when the entry is repointed at the fork. The forks inherited no upstream tags. The two newest commits on the branch are from 2026-09-07 (Xiaojie Hao); confirm that this is the intended state.

ismm: the branch leaves EC-Earth4 `main` at e95e979 (2026-06-05), after release 4.1.8 (2026-05-28) and before 4.1.9, plus our commits; hence `4.1.8-awi.1`. PISM commit confirmed as the intended state (2026-10-05).

dEBM: the `esm_tools` branch is upstream `v1.0` plus 25 commits. A later fix in a fork becomes `-awi.2`, and so on.

## FESOM: state of `main` for `2.8.1`

Checked 2026-10-05 late morning. `main` is 44e4f9c0; nothing further has to be merged before the tag.

Merged on 2026-10-05:

- **#1112** hemispheric K_GM_max. The one the tuned configuration needs.
- **#1113** Intel `-prec-sqrt` fix. Compile flags only (`src/CMakeLists.txt`); results may change at round-off level.
- **#1114** remove the DIC_PI switch and the obsolete top-level REcoM namelists. The namelists esm_tools stages (`config/bin_2p3z2d/namelist.recom` and `namelist.tra`) are untouched and never carried DIC_PI. It does move the REcoM submodule from 8f38cc4 to f71285c, so the carbon-cycle variant gets REcoM code that PICAL_cc did not run.
- **#1111** NaN-mask cavity nodes in native 2D output. CORE3 has cavities, so 2D fields are now NaN under the ice shelves where they held values before. The evaluation scripts and reval have not been run against that. (Its change to `config/bin_2p3z2d/namelist.io` does not reach us; esm_tools uses its own `namelist.io.recom`.)

Not needed:

- **#1103** IDEMIX horizontal smoothing independent of dt. Closed unmerged, not in `main`, and not in any binary of the tuning line (gmhemi_636b26aa, ccgmhemi_90074f13 and the current runs have no `idemix_hor_smooth_dtref`). Every 1800 s run with IDEMIX on, including gmhemi1800 and the IDEMIX control arm, ran without it. It is inactive without IDEMIX, and with IDEMIX it would raise the mixing above what was tuned and evaluated; the test with the fix (nx1800) did not restore Labrador convection either. So neither a TKE-only nor an IDEMIX release wants it. The long-term fix discussed in the pull request is IDEMIX2.

Still open, not for 2.8.1: #1059 (CORE3 as default mesh), #1050 (coupling interface rework), drafts #958, #897, #989, #1107, and #734, #607.

In `main` but never run on the tuning line (the production binary is an older base plus #1060, #1087 and the GM change):

- #1056 flux-form FCT low-order tendency in both precisions (changed FESOM's own reference values)
- #1109 virtual BGC tracer fluxes, #1101 threaded aEVP, #1113, #1114 with the newer REcoM, #1111, and the rest of the 140 commits since `2.8.0`

## OpenIFS and LPJ-GUESS: commits the tuning line did not have

The physics-only build tree was behind its remotes until 2026-10-05 (now fast-forwarded; its installed binaries are still the old build).

- OpenIFS, 5 commits: four reworking the surface-pressure (dry mass) fixer to compute its global norms and correction in double precision (6847539, b40698f, ddcabf0, f6d1860) and the option to send 2 m temperature to the ice sheet (d5da841). The fixer change is NOT limited to the carbon-cycle variant: it is unconditional code in `pfixer`/`qmfixer`, the fixer is on in every variant (`LMASCOR = .true.` in the physics-only fort.4 as well, esm_tools `massfixer: 2`), and OpenIFS is built single precision, which is the case the change addresses. PICAL_cc ran with it; the physics-only tuning line did not.
- LPJ-GUESS, 1 commit: soil initialisation fix for cold starts (179fb0d).

## esm_tools side

Pattern of the v3.4.2 pinning commit (ff31922d9): a version entry per component whose `branch:` is the tag, a coupling file listing the six pinned components, and the version added to `available_versions` and `choose_version` in `configs/setups/awiesm3/awiesm3.yaml`.

For v3.5.0:

- new coupling `configs/couplings/awiesm3_v3.5.0/` with `fesom-2.8.1`, `oifs-48r1v5`, `lpj_guess-4.1.13`, `oasis3mct-5.3`, `xios-2.5.4`, `rnfmap-v1.4`; the ice-sheet coupling adds the PISM, dEBM and ismm tags and repoints PISM and dEBM at the forks
- `oifs-48r1v5` currently resolves to the branch `awiesm3-develop`; change its value to the tag `48r1v5` and give develop its own entry for the moving branch
- `v3.5.0` (and, if in scope, `v3.5.0-cc`, `v3.5.0-is`) as new pinned versions; `v3.5` and `develop*` stay the moving line
- move the tuned settings from the runscripts into `awiesm3.yaml` (mixing scheme, hemispheric GM values, `albsn`, 1800 s step)
- before the esm_tools release, confirm with `git ls-remote` that each of the six references resolves to `refs/tags/...` only

## Order of work

1. Settle the tuning (IDEMIX control arm).
2. FESOM pull requests: done, #1112 is in `main`.
3. Move the tuned settings into `awiesm3.yaml`.
4. Build one tree from the exact candidate commits and run a validation against the tuning line.
5. Cut the tags on those commits.
6. Make the esm_tools pinning commit.
7. Tag esm_tools.

## Open decisions

- Decided 2026-10-05: one tag set for all three variants.
- Length of the validation run (20 years from the 2210 restart is the proposal, because #1056 and the newer OpenIFS and LPJ-GUESS commits are untested on this line).

## Validation (started 2026-10-05)

Fresh trees, `esm_master get` + `comp`, in `/albedo/scratch/user/jstreffi/model_codes_v350/` (work_user is at 978k of 1M files): `awiesm3-develop`, `awiesm3-develop-cc`, `awiesm3-develop-is`. All nine components check out at exactly the candidate commits.

- **PICAL_v350val**, physics-only: repeats PICAL_crunveg_tke_albsn082's 2190-2209 leg from the same restart with the candidate build. A pair over the same years.
- **PICAL_cc_v350val**, carbon cycle: repeats PICAL_cc's 2140-2144 from the same start. A pair over the same years; what differs is FESOM `main` 44e4f9c0 against 90074f13, and the esm_tools environment (REcoM XIOS group from the branch instead of the hand patch).
- Ice-sheet variant: built on albedo; its runscripts (`runscripts/awiesm3/develop-is/`) are levante-only, so a run needs a port or has to happen on levante.

Found on the way: esm_tools `v3.5-cc` sets `-DUSE_SINGLE_PRECISION=ON` for FESOM, while `develop-cc` sets it OFF with the comment that REcoM is not single-precision ready. To be fixed in the esm_tools release work.
