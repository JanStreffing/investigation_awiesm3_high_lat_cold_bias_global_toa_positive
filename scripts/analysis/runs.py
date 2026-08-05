"""The campaign run list, in one place.

Every evaluator used to carry its own copy of RUNS, and they drifted three times --
each time leaving one script unable to reproduce numbers the others had printed.
Add a run HERE and all evaluators pick it up.

Also holds the shared evaluation window and paths, for the same reason: Y0/Y1 were
edited from (1872,1875) to (1872,1915) by hand in three files when the runs were
extended, which is exactly the kind of edit that gets missed in one of them.

Conventions
-----------
* label      short, <=16 chars, appears in table headers
* directory  under RT; evaluators skip it with a warning if output is incomplete,
             which is the intended behaviour for in-flight or non-tuning runs
* ORDER MATTERS -- tables are printed in this order, and 'control' must be present
  because every delta is taken against it.
"""

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
OBS = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
ERA5_T2M = '/work/ab0246/a270092/obs/era5/netcdf/T2M.nc'

# 44-year window: 1870-1871 discarded as spin-up.  Detection threshold on Siberian
# JJA T2m is +-0.242 K here; it was +-0.89 K at the original 4 years.
Y0, Y1 = 1872, 1915

RUNS = [
    # ORDER IS THE ORIGINAL EVALUATOR ORDER -- do not regroup, it sets table column
    # order and keeps output comparable with tables already in the logbook/report.
    ('control', 'amip_pi_base'),
    ('A1a ovl=0.10', 'amip_A1_overlap01'),
    ('A1b ovl=0.35', 'amip_A1_overlap035'),
    ('A2 KKland=150', 'amip_A2_kknumland150'),
    ('expA rvrs=500', 'amip_expA_rvrsmin500'),
    ('A1c depth1500', 'amip_A1c_depliqdepth1500'),
    ('B1 detrpen.45', 'amip_B1_detrpen045'),
    ('B2 convi=25', 'amip_B2_clddiffconvi25'),
    ('AB ovl+convi', 'amip_AB_ovl035_convi25'),
    ('B3 clddiff', 'amip_B3_clddiff15e6'),
    ('B4 entshalp3', 'amip_B4_entshalp3'),
    ('B5 capdcycl0', 'amip_B5_capdcycl0'),
    ('B6 lcritsnow', 'amip_B6_lcritsnow1e5'),
    ('B7 rvice.22', 'amip_B7_rvice022'),
    ('B8 lamsk5', 'amip_B8_lamsk5'),
    ('ABB8 A1b+B2+B8', 'amip_ABB8'),
    ('C1 rlam75', 'amip_C1_rlam75'),
    ('C2 rlam40', 'amip_C2_rlam40'),
    ('E1 lamsk2.5', 'amip_E1_lamsk25'),
    ('D1 capdcycl4', 'amip_D1_capdcycl4'),
    ('D2a inpsea.2', 'amip_D2a_inpsea02'),
    ('D2b inp+p700', 'amip_D2b_inpsea02_p700'),
    ('piCTRL 1850', 'amip_picontrol'),
    ('F1 z0h/10', 'amip_F1_z0h10'),
    ('F2 LAI=3', 'amip_F2_lai3'),
    ('F3 cov=0.7', 'amip_F3_cov07'),
    ('F4 rsmin1000', 'amip_F4_rsmin1000'),
    ('F5 all four', 'amip_F5_allveg'),
    ('G1 F4+D2b', 'amip_G1_F4_D2b'),
    # Round 13 (H): snow cover fraction RQSNCR 1/10 -> 1/30.  REJECTED 2026-08-02 --
    # JJA +0.020 K (predicted +0.2..+0.7) and DJF -1.233 K, the worst winter damage in
    # the campaign.  Kept in RUNS because a rejected lever is still evidence.
    # Diagnosis: the albedo response lands in Sep-Nov, not June -- shallow snow is an
    # autumn phenomenon.  See monthly_lever_check.py.
    ('H1 snowcr30', 'amip_H1_snowcr30'),
    ('H2 G1+snowcr', 'amip_H2_G1_snowcr30'),
    # Round 14, complete 2026-08-03.  All namelist-only (&NAMSURFTUNE), one binary.
    # G2/G3 bracket F4: 500/1000/2000 -> +0.336/+0.521/+0.876, ACCELERATING, no knee,
    # and G3 costs DJF -0.798 -> cap RVRSMIN at 1000.  G4 adds tundra (type 9, 25.6% of
    # the box, as-released 80 s/m = lowest of any vegetated type) and is the CAMPAIGN
    # BEST: JJA +0.952 (t=7.68), SO SW RMSE 4.800, seasonally clean.  Its melt mechanism
    # is falsified though -- May/June SWE did not move.  See RUNS_AND_PARAMETERS.md.
    ('G2 rsmin500', 'amip_G2_rsmin500'),
    ('G3 rsmin2000', 'amip_G3_rsmin2000'),
    ('G4 tundra225', 'amip_G4_tundra'),
    # Round 15: sub-grid snow depletion (Niu & Yang 2007), namelist-only.
    # I1 = G4 + mode 1; I2 = mode 1 alone (isolates the snow route from the
    # vegetation route); I3 = G4 + mode 2, the SDOR-scaled scale-aware variant,
    # calibrated to match I1 at TCO95 so any difference is spatial structure.
    ('I1 scf', 'amip_I1_scf'),
    ('I2 scf only', 'amip_I2_scf_only'),
    ('I3 scf sdor', 'amip_I3_scf_sdor'),
    # Round 16: exposed-snow skin conductivity ZSNOW 7 -> 15/25, on top of I1.
    # Round 15 showed the snow-depletion fix removes a COMPENSATING error --
    # excess cover had been propping DJF up ~2.7 K, hiding a -4.7 K bias vs ERA5.
    # DJF downward LW is already correct (+0.9 W/m2 vs CERES), so the deficit is
    # non-radiative: skin decoupling in polar night.
    ('J1 lamsk15', 'amip_J1_lamsk15'),
    ('J2 lamsk25', 'amip_J2_lamsk25'),
    # Round 17: moss/litter insulating layer as a skin-conductivity proxy,
    # lambda_sk 10 -> 2.9 / 1.7 on vegetation types 4 (larch) and 9 (tundra)
    # only -- the two most boreal-confined types (62%/70% of their global area
    # is >50N).  Gaillard 2025 attributes most of their >2 K high-latitude
    # summer gain to a surface organic layer, which HTESSEL lacks entirely.
    ('L1 moss2.9', 'amip_L1_moss29'),
    ('L2 moss1.7', 'amip_L2_moss17'),
    # Round 17 (K): the snow-free LAND ALBEDO residual, on the G4 base.  Global
    # land albedo is +0.0154 too high vs CERES; masking snow per cell per month
    # splits that into +0.0074 snow (the I-series lever) and +0.0080 snow-free.
    # Every FOREST type is already correct once snow is masked (all within
    # +0.006), so RVVEGALB's high-vegetation entries are untouched; the error is
    # on sparse/cultivated surfaces.  K1 corrects the vegetation table (crops
    # x0.887, irrigated crops x0.910, semidesert x0.906) AND the bare-soil
    # background (x0.95); K2 is the soil background ALONE, isolating the ~25% of
    # land that RVVEGALB cannot reach (its table row is identically zero there).
    # TUNDRA IS DELIBERATELY EXCLUDED so K stays orthogonal to J and L, which both
    # act on boreal tiles.  NOT a Siberian lever -- score on bias_by_tile.py and
    # global T2m RMSE, with the Siberian box as a guardrail only.
    ('K1 landalb', 'amip_K1_landalb'),
    ('K2 soilalb', 'amip_K2_soilalb'),
]

# Not tuning levers -- do not add these to RUNS.  Eleven LPJG forcing-generator runs
# from 3-5 July 2026 emit the daily LPJG forcing set rather than evaluation fields, so
# they legitimately have no atm_remapped_1m_2t_* output, plus one superseded A2 attempt.
# See the run inventory in notes/ATMOSPHERE_TUNING_LOGBOOK.md.
NOT_LEVERS = [
    'amip_A2_kkland150',
    'amip_lpjgforce_chk', 'amip_nolpjg_forc', 'amip_nolpjg_forcing', 'amip_nolpjg_pi1870',
    'amip_pi_clean1', 'amip_pi_clean2', 'amip_pi_dbg1', 'amip_pi_dbg2', 'amip_pi_dbg3',
    'amip_pi_fixtest', 'amip_pi_forcing',
    # Round 18 aerosol DIAGNOSTICS (2026-08-05). These sit on the PRESENT-DAY base
    # (1989-2015), not amip_pi_base (1872-1915), because MACv2-SP is transient and
    # at 1870-1915 the anthropogenic plumes are already near zero -- the test would
    # show nothing there. They must NOT enter the scored RUNS list: evaluate.sh
    # would difference them against the PI control over the wrong years.
    # Compare them against amip_presentday, and on BOTH reflection columns
    # (clear-sky and cloud) by latitude band, never the all-sky total.
    'amip_M1_noanthaer',   # LMACV2SP=.false. -- anthropogenic aerosol removed
    'amip_M2_aer3d',       # LAER3D=.true.    -- 3D CAMS vertical distribution
    # Round 19 N series: DAILY snow/soil process diagnostic on the K1 base, 10 yr only
    # (1870-1880). Not tuning candidates and far too short for the 44-yr thresholds --
    # they exist to resolve whether the winter soil collapse is seeded by the October
    # cover deficit, and whether the pack is ripe when the spring depletion is needed.
    # N2 minus N1 isolates ECE_SNOW_SCF=1 at daily resolution.
    'amip_N1_snowdiag',
    'amip_N2_snowdiag_scf',
    # Round 19 O series: the SAME scheme re-parameterised so density actually separates
    # autumn from spring (z0 0.018, rho_new 170, m 4.0 / 3.0 vs the current 0.016/100/1.6).
    # Namelist-only.  Same K1 base, 10 yr and daily output as N1/N2 so all four difference
    # directly.  Not scorable against the 44-yr thresholds -- diagnostic pair, not levers.
    'amip_O1_scf_m4',
    'amip_O2_scf_m3',
]
