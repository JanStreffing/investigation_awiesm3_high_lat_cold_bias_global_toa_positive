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
    # Round 19: the snow-cover depletion revisited on the K1 base, at FULL campaign
    # length so it is scorable against the 44-yr thresholds.  N2 is the scheme at its
    # current parameters (z0 0.016, rho_new 100, m 1.6); O1/O2 re-parameterise it so
    # density actually separates autumn from spring.  The ratio that matters is
    #     x_Oct/x_May = (d_Oct/d_May)*(rho_May/rho_Oct)^m = 0.469*1.853^m
    # so m>1.23 gives the right ordering but m=1.6 separates the two regimes by only
    # 0.077 -- which is why the -0.075 October cover deficit that seeds the winter soil
    # collapse buys a spring depletion just 0.08 deeper.  m=4 gives a factor 5.5, but
    # rho_new must be recentred on autumn density first or October depletes harder still.
    # All three carry daily sd/rsn/tsn/asn/stl1/stl2 via a per-run file_def override.
    ('N2 scf current', 'amip_N2_snowdiag_scf'),
    ('O1 scf m4', 'amip_O1_scf_m4'),
    ('O2 scf m3', 'amip_O2_scf_m3'),
    # Round 20: the depletion curve REBUILT on observations after the tanh was
    # falsified outright.  tanh(x) has range (0,1) OPEN AT 1, so complete snow cover
    # is not representable at any (z0, rho_new, m) -- while RIHMI-WDC snow courses
    # report exactly 10/10 in 99.74% of 14369 Siberian DJF surveys, and the tanh
    # reaches complete cover on 29% of January field cases against an observed 99.6%.
    # Cost, measured against 174676 station observations at 0.2 m (-6.1 degC):
    # N1 (off) -5.3 = +0.8 bias, N2 -26.3 = -20.2, O1 -16.1, O2 -14.7.  ECE_SNOW_SCF=3
    # replaces it with min(1,(d/d_c)^b), d_c = SCALE*DC*(rho/200)^4.7, split by
    # vegetation type and fraction, fitted to 36492 snow-course surveys.  P1 takes the
    # fit literally (SCALE=1); P2 scales d_c x3 for 100 km sub-grid variance, the one
    # parameter station data cannot constrain.  See notes/RUNS_AND_PARAMETERS.md.
    # Round 20b: P1/P2 are WITHDRAWN and replaced by P3/P4.  P1 aborted at 1888-03-22
    # with ABOR1 'Very snow cold temperature' (srfsn_webal_mod.F90:451, PTSN < 100 K):
    # Tsn 91.9 K, SWE 0.056 kg/m2, PFRSN 1.000.  The fitted power law was extrapolated
    # below any depth the snow courses sampled -- at rho=100 d_c = 5.4e-4 m, so a
    # 0.56 mm dusting got FULL COVER and the snow tile took the whole grid-box flux
    # into ~zero heat capacity (SCF/SWE 17.9 vs 0.10 as-released).  ECE_SNOW_SCF_SWEMIN
    # now floors d_c at SWEMIN/rho -- a minimum snow MASS, the right currency for a
    # heat-capacity failure.  P2 never aborted but carries the same defect, so its
    # output is tainted too and it was cancelled at leg 5.
    ('P3 scf fit', 'amip_P3_scffit'),
    ('P4 scf fit x3', 'amip_P4_scffit_x3'),
    # Round 21: SWEMIN calibrated against Rutgers.  P3 cured the winter but took
    # Siberian September cover to +0.259 against the satellite record, worse than the
    # as-released ramp's +0.107, because the fitted d_c for fresh low-density snow is
    # ~2 mm and a dusting saturates.  SWEMIN floors d_c at SWEMIN/rho and binds only at
    # low density, so it is an almost pure autumn control.  An offline sweep scored over
    # the whole Sep-May season put the optimum at 30 (RMSE 0.045 vs P3's 0.139) with an
    # interior maximum -- past 30 the September gain is paid for by Nov/Dec going
    # deficient.  P5/P6 bracket the albedo->melt->SWE feedback the sweep holds fixed.
    # NOTE both ran on the pre-DCMAX-fix binary (562df81), where the 0.30 m cap also
    # clipped the floor; reconstructions must use the old operator order.
    ('P5 swemin15', 'amip_P5_swemin15'),
    ('P6 swemin30', 'amip_P6_swemin30'),
    # Round 23 (S): the Southern Ocean cloud-AMOUNT deficit, after round 22 showed it is
    # 65.5% of a +7.46 W/m2 band error and unreachable through the CRE/opacity route.
    # S1/S2/S3 were WITHDRAWN before completion (see NOT_LEVERS); S4 is the survivor.
    # RCL_INPPMIN 70000 -> 50000 extends D2b's marine ice-nuclei floor from below ~700 hPa
    # up to ~500 hPa, i.e. more of the mixed-phase column, while still sparing tropical
    # anvils -- which is the whole point, because every global low-cloud knob overshoots a
    # tropical band that is only -0.67 W/m2 from CERES period-clean.
    # ADDED 2026-08-10: the run completed 48 years on 2026-08-09 and had never been
    # entered here, so evaluate.sh had never scored it.
    ('S4 inppmin50k', 'amip_S4_inppmin50000'),
    # Round 27 (2026-08-10), the L-series: the first full-length runs to carry RSNOWLIN2.
    #
    # RSNOWLIN2 0.030 -> 0.04 is the temperature exponent of ice->snow autoconversion
    # (ZZCO=PTSPHY*RSNOWLIN1*EXP(RSNOWLIN2*(T-RTT)), cloudsc.F90:2858).  Raising it slows
    # the conversion 31 % in 200 K cirrus but only 11 % at 250 K, so it builds cold thin
    # cirrus rather than thickening anvils -- and thin cirrus is LONGWAVE-dominant.  At one
    # year it gave +2.74 W/m2 of tropical LW CRE, 124 % of the -2.16 deficit and five times
    # the largest response among the other 50 arms, with tropical and global net TOA both
    # below their thresholds.  The value 0.04 is the awiesm3 v3.4.2 TCO95L91-DARS2 setting.
    #
    # THE CONTROL FOR ALL FOUR IS P5, not amip_pi_base.  Verified 2026-08-10 that T3, the
    # base these were built on, is BIT-FOR-BIT identical to P5 at year 1870 (SO SW CRE,
    # tropical LW CRE, global net TOA all to 4 dp); the only namelist difference is
    # ECE_CLIMR_DMS=.true., which reads the DMS field but does nothing at S=0.
    #
    # WHAT THE ONE-YEAR SCREENS COULD NOT ANSWER, and why these exist: Siberia.  The JJA
    # threshold is +-0.242 K and DJF +-0.588 K, both 44-year numbers.  Every arm here
    # changes cloud, and LY2 uses RCL_OVERLAPLIQICE=0.1 -- A1a's value, rejected in round 3
    # for boreal damage.  These runs are here to price that.
    ('LX1 rsnow', 'amip_LX1_long'),           # RSNOWLIN2 alone; the gating run
    ('LX3 rsnow+dms+inp', 'amip_LX3_long'),   # + DMS S=166 + RCL_INPPMIN 50000
    ('LY2 rsnow+ovl0.1', 'amip_LY2_long'),    # + RCL_OVERLAPLIQICE 0.10

    # Round 30 (2026-08-19).  All three are S4 + ONE namelist number, 46 yr, so S4 is
    # the control and no new control was needed.  They exist to answer two questions the
    # campaign had assumed rather than measured.
    #
    # N1/N2 -- has RCL_INPSEA any headroom?  It has sat at 0.2 in all 46 runscripts that
    # set it and was NEVER scanned; 0.2 is justified as the "marine biogenic floor", an
    # argument about INP concentration.  Measured against cloud phase instead
    # (so_cloud_phase_vs_goccp.py): the model's low-level SO liquid fraction is 0.875
    # against CALIPSO-GOCCP's 0.858, i.e. already MORE liquid than observed.  So these
    # are expected to move SO CRE and to be REJECTED anyway, for pushing phase past the
    # observation.  Scored to confirm that rather than assert it.
    #
    # W1 -- is round 23's tropical cost real?  This IS S3, designed then CANCELLED on an
    # inference by analogy with B3/RCLDIFF, a GLOBAL cloud-erosion term.  RCLCRIT_SEA acts
    # on stratiform warm-rain autoconversion over sea (PLSM > 0.5), which is the Southern
    # Ocean's regime while tropical precipitation is largely convective.  SCORE ON CLOUD
    # AREA, not CRE -- round 23's instruction, because CRE conflates the amount and
    # opacity thirds and the INP branch already looked good on CRE while doing nothing
    # for area.  Disqualifier: tropics beyond +-0.5.
    ('N1 inpsea0.1', 'N1'),                   # RCL_INPSEA 0.2 -> 0.10
    ('N2 inpsea0.05', 'N2'),                  # RCL_INPSEA 0.2 -> 0.05
    ('W1 clcritsea6e4', 'W1'),                # RCLCRIT_SEA 2.5e-4 -> 6.0e-4
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
    # N1 is the daily-output TWIN of K1 (identical configuration, scheme off).  Kept out
    # of RUNS so it does not appear as a duplicate column; it doubles as a reproducibility
    # check -- N1 must match K1 to within the noise floor.
    'amip_N1_snowdiag',
    # Round 20 first attempt, WITHDRAWN 2026-08-06 -- unfloored d_c gave full snow
    # cover to a sub-millimetre pack.  P1 aborted at 1888-03-22; P2 was cancelled at
    # leg 5 carrying the same defect.  Superseded by P3/P4 with ECE_SNOW_SCF_SWEMIN.
    # Kept out of RUNS so evaluate.sh cannot score partial, defective output.
    'amip_P1_scffit',
    'amip_P2_scffit_x3',
    # Round 23 (S), WITHDRAWN 2026-08-09 before completion.  S2 is a DUPLICATE of A1b
    # (RCL_OVERLAPLIQICE=0.35, already run 44 years); S1 and S3 were cancelled with it
    # when the round was re-scoped.  Only S4 survived and is in RUNS above.
    'amip_S1_clddiff3e6', 'amip_S2_ovlliqice035', 'amip_S3_clcritsea6e4',
    # Rounds 24-26, the DMS/CCN series.  ONE-YEAR screening runs, and their control is
    # T3 (S=0, same DMS input, same winds), NOT amip_pi_base.  Putting them in RUNS would
    # make evaluate.sh difference a single year against a 44-year PI control -- the exact
    # window mismatch NOT_LEVERS exists to prevent.  Scored instead by
    # dms_ccn_bracket.py, which computes its own per-band 1-year thresholds.
    'amip_T1_dmsread',      # T1: does the eighth ICMCL field read at all
    'amip_T2_dmsflux',      # T2: sea-air flux diagnostic, shipped DMS climatology
    'amip_T3_dmsrev3',      # T3: DMS-Rev3 ingested; the CONTROL for the U series
    'amip_U1_dmsccn166',    # U1: ECE_DMS_CCN_SENS = 166, top of Woodhouse's range
    'amip_U2_dmsccn43',     # U2: 43, the published floor
    'amip_U3_dmsccn90',     # U3: 90, the midpoint
    # Rounds 26-27 one-year screens.  Same reason as the T/U series: ONE year, and their
    # controls are T3 or X1, not amip_pi_base, so evaluate.sh would difference a single
    # year against a 44-year PI control.  Scored by their own scripts instead --
    # w_series_screen.py and x_series_ladder.py, which carry 1-year thresholds.
    'amip_V1_entrorg125', 'amip_V2_detrpen111',
    'amip_W1_rprcon148', 'amip_W2_entrorg207', 'amip_W3_rmfdeps048', 'amip_W4_detrpen132',
    'amip_W5_entrdd108', 'amip_W6_rvice018', 'amip_W7_lcritsnow146', 'amip_W8_rsnowlin204',
    'amip_W9_dars2stack',
    'amip_X1_stack_rsnow', 'amip_X2_stack_rsnow_dms166',
    'amip_X3_stack_rsnow_dms166_inp50k', 'amip_X4_stack_rsnow_inp50k',
    'amip_Y1_ovl035', 'amip_Y2_ovl01',
    # Z1 is a FORCING mechanism test, not a lever: it checks whether NCMIPFIXYR=1850 can be
    # pinned from the live NAMECECMIP block.  If it works it changes the absolute energy
    # baseline by ~0.15 W/m2 and is therefore NOT comparable to this archive at all.
    'amip_Z1_fixyr1850',
    # LX4 CANCELLED 2026-08-10 after ~12 min, before any output.  Dominated: it warms the
    # Southern Ocean in net (+0.377) and does nothing for the energy target (+0.067), and
    # its unique content -- RCL_INPPMIN 50000 at full length on the stack -- already exists
    # as S4 minus P5.  It was LX3 without the only lever that cools the SO.
    'amip_LX4_long',
]
