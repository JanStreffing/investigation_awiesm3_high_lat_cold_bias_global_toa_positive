"""Is (depth, density) SUFFICIENT, or does the depletion curve need a melt variable?

WHY THIS EXISTS.  A snow-cover scheme in a climate model must be a function of model
STATE, never of the calendar -- otherwise it is pinned to the present-day seasonal
cycle and cannot be trusted in a 10 K warmer world or a glacial.  The fitted curve
    SCF = min(1, (d/d_c)**b),   d_c = d_c0*(rho/100)**m
is calendar-blind by construction: its only arguments are snow depth and bulk
density.  But "calendar-blind" is not the same as "sufficient".  If the real cover
at a GIVEN (d, rho) is systematically lower in spring than in autumn, then density
does not capture the melt state, the seasonal signal in the fit is coming from the
seasonal cycle of the inputs rather than from physics, and the curve will mis-predict
in any climate whose d-rho trajectory differs from today's.

THE TEST.  Match surveys on (depth, density) and compare ACCUMULATION season
(Oct-Jan) against ABLATION season (Mar-May) within the same bin.  Under the
hypothesis that (d, rho) is sufficient, the two means agree inside noise.  A
systematic ablation deficit is direct evidence of hysteresis that density misses.

This is the same question the melt-state gate was meant to answer in round 19, but
that attempt was rejected on model output (May ripeness 0.04 -- the gate would have
fired in June).  Here it is asked of observations, where it can actually be settled.

IF HYSTERESIS IS REAL, the fix must still be a state variable, not a date.  The model
already carries two candidates as snow prognostics -- liquid water content PWSNM1M and
snow temperature PTSNM1M -- neither of which currently reaches SURFBC_CTL, so using
one means extending the argument list, not adding a prognostic.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

OBSD = '/work/ab0246/a270092/obs/RIHMI-WDC/data'
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
TYPES = {0: 'field (pole)', 1: 'forest (les)'}
ACC, ABL = (9, 10, 11, 0), (2, 3, 4)          # Oct-Jan accumulation, Mar-May ablation
DBINS = np.array([.05, .08, .12, .17, .23, .30, .40, .55, .80])
RBINS = np.array([120, 150, 180, 210, 250, 300, 400])

with xr.open_dataset(f'{OBSD}/RIHMI-WDC_snmar.nc') as ds:
    la, lo = ds['lat'].values, ds['lon'].values
    k = (la >= 55) & (la <= 75) & (lo >= 60) & (lo <= 180)
    cov = ds['fraction_of_the_snow_course_covered_by_snow'].isel(station=k).values / 10.0
    dep = ds['snow_depth_mean'].isel(station=k).values / 100.0
    rho = ds['snow_density'].isel(station=k).values * 1000.0
    tm = ds['time'].values
mo = tm.astype('datetime64[M]').astype(int) % 12

print(__doc__.split('IF HYSTERESIS IS REAL')[0])
print('=' * 96)

for t in (0, 1):
    c, d, r = cov[:, :, t], dep[:, :, t], rho[:, :, t]
    m3 = np.repeat(mo[None, :], c.shape[0], axis=0)
    ok = np.isfinite(c) & np.isfinite(d) & np.isfinite(r) & (d > 0) & (r > 30) & (r < 600)
    c, d, r, mm = c[ok], d[ok], r[ok], m3[ok]
    acc, abl = np.isin(mm, ACC), np.isin(mm, ABL)
    print(f'\n\n### TYPE {t+1} -- {TYPES[t]}   accumulation n={acc.sum()}  ablation n={abl.sum()}\n')
    print('  Mean observed cover at MATCHED (depth, density).  acc = Oct-Jan, abl = Mar-May.\n')
    print(f'  {"depth m":>12s}{"rho":>10s}{"n acc":>8s}{"n abl":>8s}{"acc":>8s}{"abl":>8s}{"abl-acc":>9s}')
    diffs, wts = [], []
    for i in range(len(DBINS) - 1):
        for j in range(len(RBINS) - 1):
            s = (d >= DBINS[i]) & (d < DBINS[i + 1]) & (r >= RBINS[j]) & (r < RBINS[j + 1])
            na, nb = (s & acc).sum(), (s & abl).sum()
            if na < 25 or nb < 25:
                continue
            ca, cb = c[s & acc].mean(), c[s & abl].mean()
            diffs.append(cb - ca); wts.append(min(na, nb))
            print(f'  {DBINS[i]:.2f}-{DBINS[i+1]:.2f}'.rjust(14)
                  + f'{RBINS[j]:>4d}-{RBINS[j+1]:<5d}'.rjust(10)
                  + f'{na:8d}{nb:8d}{ca:8.3f}{cb:8.3f}{cb-ca:+9.3f}')
    if diffs:
        diffs, wts = np.asarray(diffs), np.asarray(wts, float)
        wm = np.sum(diffs * wts) / wts.sum()
        # paired sign test across bins -- distribution free, no normality assumed
        npos, nneg = (diffs > 0).sum(), (diffs < 0).sum()
        print(f'\n  {len(diffs)} matched bins.  weighted mean (abl - acc) = {wm:+.4f}')
        print(f'  bins where ablation is LOWER: {nneg}/{len(diffs)}   higher: {npos}/{len(diffs)}')
        print(f'  unweighted mean {diffs.mean():+.4f}   median {np.median(diffs):+.4f}'
              f'   min {diffs.min():+.4f}  max {diffs.max():+.4f}')
    else:
        print('\n  no bins with enough surveys in both seasons')

print("""

  READING IT.  A weighted mean near zero with the sign split roughly even means
  (depth, density) already carries the melt state: the curve can be a pure function of
  the two variables SURFBC_CTL already receives, and the seasonal cycle emerges from
  the snow physics rather than being imposed.  A consistently NEGATIVE (abl - acc)
  across most bins means real hysteresis remains after density is accounted for, and
  the curve needs a melt-state variable -- liquid water content PWSNM1M is the
  physically correct one and is already prognostic, but it is not currently passed to
  SURFBC_CTL, so it would mean extending the argument list through the caller chain.
  Either way the answer is a state variable; the calendar is not an option.""")
