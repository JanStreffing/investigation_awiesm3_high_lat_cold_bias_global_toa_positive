"""Paired difference between two coupled arms that share a branch point.

Branched arms start from one state and diverge, so the difference series begins at zero and
grows; the unpaired 1.96*sd*sqrt(2/n) threshold discards that pairing and has hidden a
resolved result in this campaign before.  Everything here is a paired t-test on the annual
series, which is what the report's round-29 note requires.

Usage:  A=PICAL_momixoff B=PICAL_v35def Y0=2021 Y1=2049 python3 scripts/analysis/paired_arm_diff.py
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, glob, warnings
from scipy import stats
warnings.filterwarnings('ignore')

R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
A, B = os.environ.get('A', 'PICAL_momixoff'), os.environ.get('B', 'PICAL_v35def')
Y0, Y1 = int(os.environ.get('Y0', 2021)), int(os.environ.get('Y1', 2049))
SEAS = {'ANN': list(range(12)), 'DJF': [11, 0, 1], 'JJA': [5, 6, 7]}
BOXES = [('GLOBAL', -90, 90, 0, 360, None), ('60-90N', 60, 90, 0, 360, None),
         ('Siberia 55-75N', 55, 75, 60, 180, True), ('30-60N', 30, 60, 0, 360, None),
         ('20S-20N', -20, 20, 0, 360, None), ('60-90S', -90, -60, 0, 360, None)]

f0 = sorted(glob.glob(f'{R}/{A}/outdata/oifs/atm_remapped_1m_lsm_*.nc'))[0]
with xr.open_dataset(f0, decode_times=False) as d:
    msk = np.squeeze(d['lsm'].values); msk = msk[0] if msk.ndim == 3 else msk
    lat = np.squeeze(d['lat'].values); lon = np.squeeze(d['lon'].values)
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], msk.shape)
LAT = np.broadcast_to(lat[:, None], msk.shape); LON = np.broadcast_to(lon[None, :], msk.shape)


def series(arm):
    """annual (year, month, lat, lon) T2m"""
    out = []
    for y in range(Y0, Y1 + 1):
        g = (glob.glob(f'{R}/{arm}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc')
             or glob.glob(f'{R}/{arm}/outdata/oifs/atm_remapped_1m_2t_1m_{y}-{y}.nc'))
        if not g:
            raise SystemExit(f'{arm}: missing year {y}')
        with xr.open_dataset(g[0], decode_times=False) as d:
            out.append(np.squeeze(d['2t'].values))
    return np.stack(out)


SA, SB = series(A), series(B)
print(__doc__.split('Usage:')[0])
print(f'{B} minus {A}, {Y0}-{Y1} (n={Y1-Y0+1}), paired by year\n')
print(f'{"box":<18}{"season":>7}{A[:13]:>14}{B[:13]:>14}{"diff":>9}{"t":>7}{"p":>9}')
for name, la, lb, lo, hi, land in BOXES:
    # hi-lo == 360 is the whole globe; (hi-lo) % 360 would be 0 and select one meridian.
    k = (LAT >= la) & (LAT <= lb)
    if (hi - lo) % 360 != 0 or hi == lo:
        k = k & (((LON - lo) % 360) <= ((hi - lo) % 360))
    if land:
        k = k & (msk > 0.5)
    for s in ('ANN', 'JJA', 'DJF'):
        a = np.array([float(np.average(SA[i][SEAS[s]].mean(0)[k], weights=W[k])) for i in range(SA.shape[0])])
        b = np.array([float(np.average(SB[i][SEAS[s]].mean(0)[k], weights=W[k])) for i in range(SB.shape[0])])
        t, p = stats.ttest_rel(b, a)
        star = '*' if p < 0.05 else ' '
        print(f'{name:<18}{s:>7}{a.mean()-273.15:14.3f}{b.mean()-273.15:14.3f}'
              f'{b.mean()-a.mean():+9.3f}{t:7.2f}{p:9.4f}{star}')
print('\n  * p<0.05 on the paired test.  Values in degC; diff in K.')
