"""Year-by-year regional state of a coupled arm: is a band settled, or still moving?

WHY THIS EXISTS.  A 20-year mean cannot tell an equilibrated bias from the middle of a
drift, and this campaign has twice adopted a number that was a trajectory.  The albedo
runaway was invisible in extent and only obvious in a per-year volume series.  So every
band here is printed per year AND as a least-squares trend with a 95 % interval, so a
"bias" and a "drift" cannot be confused.

TRAPS (campaign protocol, see the awiesm3 skill).
  * IFS TOA fluxes are ACCUMULATED J/m2 per output step -> divide by 3600.  Guarded by
    asserting incoming solar lands near 340 W/m2.
  * Sea ice EXTENT (concentration > 0.15) is not AREA.  Extent only here, labelled.
  * ERA5/CERES are present-day; a PI arm is supposed to sit below them.  This script
    prints the model's own numbers and TRENDS, not biases, so the epoch offset does not
    enter.  Use pi_offset_vs_era5.py for the bias question.

Usage:  ARMS=PICAL_ccnice,PICAL_momixoff Y0=1941 Y1=1969 \
          python3 scripts/analysis/ccnice_regional_drift.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R    = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
LSMF = ('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/'
        'outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
ACC  = 3600.0
Y0, Y1 = int(os.environ.get('Y0', 1941)), int(os.environ.get('Y1', 1969))
ARMS = [a.strip() for a in os.environ.get('ARMS', 'PICAL_ccnice').split(',')]
JJA  = [5, 6, 7]

with xr.open_dataset(LSMF, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values); m = m[0] if m.ndim == 3 else m
    lat = np.squeeze(d['lat'].values); lon = np.squeeze(d['lon'].values)
land, ocean = m > 0.5, m <= 0.5
W   = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], m.shape).copy()
bnd = lambda lo, hi: np.broadcast_to(((lat >= lo) & (lat < hi))[:, None], m.shape)
SIB = bnd(55, 75) & np.broadcast_to(((lon >= 60) & (lon <= 180))[None, :], m.shape) & land
CELL = np.broadcast_to((6.371e6**2 * np.cos(np.deg2rad(lat))
                        * np.deg2rad(abs(lat[1]-lat[0])) * 2*np.pi/m.shape[1])[:, None],
                       m.shape)

def am(f, s):
    k = s & np.isfinite(f)
    return float(np.average(f[k], weights=W[k])) if k.any() else np.nan

def yr(arm, var, y):
    p = f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p, decode_times=False) as d:
        k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values)               # (12, lat, lon)

def trend(years, vals):
    """least squares slope per decade with a 95 % interval from the residual scatter"""
    x = np.asarray(years, float); y = np.asarray(vals, float)
    k = np.isfinite(y)
    if k.sum() < 4: return np.nan, np.nan
    x, y = x[k], y[k]; n = len(x)
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b*x)
    se = np.sqrt((resid**2).sum() / (n-2) / ((x - x.mean())**2).sum())
    return b*10.0, 1.96*se*10.0

METRICS = ['net TOA', 'T2m global', 'T2m 60-90N', 'T2m 60-90S', 'T2m JJA Siberia',
           'SW CRE 45-65S', 'NH ext Mar', 'NH ext Sep', 'SH ext Feb', 'SH ext Sep']

print(__doc__.split('\n')[0]); print(f'window {Y0}-{Y1}, arms {ARMS}\n')
for arm in ARMS:
    series = {k: [] for k in METRICS}; years = []
    for y in range(Y0, Y1+1):
        t2 = yr(arm, '2t', y)
        if t2 is None: continue
        t2 = t2 - 273.15 if np.nanmean(t2) > 100 else t2
        tsr, ttr, tsrc, ci = (yr(arm, v, y) for v in ('tsr', 'ttr', 'tsrc', 'ci'))
        years.append(y)
        series['net TOA'].append(am((tsr+ttr).mean(0)/ACC, np.ones_like(m, bool))
                                 if tsr is not None else np.nan)
        series['T2m global'].append(am(t2.mean(0), np.ones_like(m, bool)))
        series['T2m 60-90N'].append(am(t2.mean(0), bnd(60, 90)))
        series['T2m 60-90S'].append(am(t2.mean(0), bnd(-90, -60)))
        series['T2m JJA Siberia'].append(am(t2[JJA].mean(0), SIB))
        series['SW CRE 45-65S'].append(am((tsr-tsrc).mean(0)/ACC, bnd(-65, -45))
                                       if tsrc is not None else np.nan)
        if ci is not None:
            nh, sh = bnd(0, 90) & ocean, bnd(-90, 0) & ocean
            e = lambda mo, msk: float(np.nansum(((ci[mo] > 0.15)*CELL)[msk]))/1e12
            series['NH ext Mar'].append(e(2, nh)); series['NH ext Sep'].append(e(8, nh))
            series['SH ext Feb'].append(e(1, sh)); series['SH ext Sep'].append(e(8, sh))
        else:
            for k in ('NH ext Mar','NH ext Sep','SH ext Feb','SH ext Sep'):
                series[k].append(np.nan)

    if not years: print(f'{arm}: no complete years in window\n'); continue
    print(f'=== {arm}  ({years[0]}-{years[-1]}, {len(years)} yr) ===')
    print(f'{"metric":<18}{"first5":>9}{"last5":>9}{"delta":>9}'
          f'{"trend/decade":>15}{"95% CI":>10}   verdict')
    for k in METRICS:
        v = np.asarray(series[k], float)
        if not np.isfinite(v).any(): continue
        f5, l5 = np.nanmean(v[:5]), np.nanmean(v[-5:])
        b, ci95 = trend(years, v)
        sig = 'DRIFTING' if np.isfinite(ci95) and abs(b) > ci95 else 'settled'
        print(f'{k:<18}{f5:>9.3f}{l5:>9.3f}{l5-f5:>9.3f}{b:>15.3f}{ci95:>10.3f}   {sig}')
    print()
