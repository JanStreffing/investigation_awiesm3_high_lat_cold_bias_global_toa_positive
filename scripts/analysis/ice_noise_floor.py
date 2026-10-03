"""Interannual noise on the sea-ice volume metrics, and what a 10-year leg can resolve.

A trend is only evidence above the scatter it is fitted through, and this campaign has
twice adopted a number that was noise.  This detrends a quiet stretch, takes the residual
sd as the 1-year noise, and converts it into the smallest trend a 10-year leg could
distinguish from zero.

For a least-squares slope over n evenly spaced years with residual sd sigma,
se(slope) = sigma / sqrt(sum (x-xbar)^2); for n = 10 that sum is 82.5, so the 95 %
detection threshold on a per-decade trend is 1.96 * 10 * sigma / 9.083 = 2.16 * sigma.

Usage:  ARM=PICAL_ccnice Y0=1941 Y1=1969 python3 scripts/analysis/ice_noise_floor.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R   = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PICAL_ccnice')
Y0, Y1 = int(os.environ.get('Y0', 1941)), int(os.environ.get('Y1', 1969))
N = int(os.environ.get('NLEG', 10))

def v(var, y):
    p = f'{R}/{ARM}/outdata/fesom/{var}.fesom.{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p, decode_times=False) as d:
        a = np.squeeze(np.asarray(d[var].values, float))/1e3   # units 1e9 m3 -> 1e3 km3
    return a if a.size == 12 else None

yrs, NH, SH = [], [], []
for y in range(Y0, Y1+1):
    n, s = v('sivoln', y), v('sivols', y)
    if n is None or s is None: continue
    yrs.append(y); NH.append(n); SH.append(s)
NH, SH, x = np.array(NH), np.array(SH), np.array(yrs, float)

rows = [('NH annual', NH.mean(1)), ('NH Apr', NH[:, 3]), ('NH Sep', NH[:, 8]),
        ('SH annual', SH.mean(1)), ('SH Sep', SH[:, 8]), ('SH Feb', SH[:, 1])]

xx = np.arange(N, dtype=float)
denom = ((xx - xx.mean())**2).sum()          # 82.5 for N = 10
fac = 1.96 * 10.0 / np.sqrt(denom)

print(__doc__.split('\n')[0]); print(f'{ARM} {yrs[0]}-{yrs[-1]}, {len(yrs)} yr, leg length {N}\n')
print(f'{"metric":<12}{"sigma (detrended)":>19}{"1-yr 95% diff":>15}'
      f'{"min trend/dec":>15}   [10^3 km3]')
for name, ser in rows:
    b, a = np.polyfit(x, ser, 1)
    sig = (ser - (a + b*x)).std(ddof=2)
    print(f'{name:<12}{sig:>19.3f}{1.96*np.sqrt(2)*sig:>15.3f}{fac*sig:>15.3f}')
print('\n  sigma is the residual scatter after removing the linear trend of the window.')
print('  "min trend/dec" is the smallest per-decade trend a 10-year leg resolves at 95 %.')
