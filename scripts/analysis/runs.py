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
]
