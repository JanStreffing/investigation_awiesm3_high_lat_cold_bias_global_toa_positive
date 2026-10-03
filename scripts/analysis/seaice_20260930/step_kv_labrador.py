"""How much weaker is vertical mixing at 1800 s in the Labrador Sea, month by month?

Globally the first-month Kv/Av below 100 m is ~4 % lower at 1800 s than at 1200 s, on both
builds (step_diff_first_month.txt).  Here: the same comparison inside the Labrador interior
(56-62N, 60-50W), per month of the first branch year, as the ratio of area-weighted mean Kv
(1800 s / 1200 s) in depth bands.  Pairs share start state and build.  The noise pair
(same build, 1 vs 4 threads) gives the reference for chance.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4'
PAIRS = [('noise 1v4 thr', ('TEST_nolockB1', 2100), ('TEST_nolockB4', 2100)),
         ('perfE/perfF old', ('TEST_perfE', 2100), ('TEST_perfF', 2100)),
         ('ob1200/ob1800 old', ('PICAL_crunveg_ob1200', 2120), ('PICAL_crunveg_ob1800', 2120)),
         ('ts1200/ihf0 new', ('PICAL_crunveg_ts1200', 2120), ('PICAL_crunveg_ihf0', 2120))]
BANDS = ((0, 50), (50, 200), (200, 600), (600, 1500))
MN = 'JFMAMJJASOND'


def load(exp, y, v):
    with xr.open_dataset(f'{R}/{exp}/outdata/fesom/{v}.fesom.{y}.nc', decode_times=False) as d:
        a = d[v].values.astype('f8'); dz = [c for c in d[v].dims if c.startswith('nz')][0]
        dep = np.abs(d[dz].values); lo = ((d['lon'].values + 180) % 360) - 180; la = d['lat'].values
    k = (la >= 56) & (la < 62) & (lo >= -60) & (lo < -50)
    if a.shape[1] == len(dep): a = np.swapaxes(a, 1, 2)          # (t, node, depth)
    a = a[:, k, :]
    return np.where(np.abs(a) > 1e30, np.nan, a), dep


for lab, (ea, ya), (eb, yb) in PAIRS:
    for v in ('Kv', 'Av'):
        a, dep = load(ea, ya, v); b, _ = load(eb, yb, v)
        nt = min(a.shape[0], b.shape[0])
        print(f'{lab:18s} {v}  ratio 1800/1200 (noise: 4thr/1thr) of Labrador-interior mean, by month')
        for z0, z1 in BANDS:
            kz = (dep >= z0) & (dep < z1)
            r = [np.nanmean(b[m][:, kz]) / np.nanmean(a[m][:, kz]) for m in range(nt)]
            print(f'   {z0:4d}-{z1:4d} m  ' + '  '.join(f'{MN[m]}:{x:5.2f}' for m, x in enumerate(r)), flush=True)
