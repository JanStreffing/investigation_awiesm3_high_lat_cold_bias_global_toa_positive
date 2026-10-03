"""Score the 2-month time-step screens: which change removes the 1800 s mixing deficit?

Metric: global mean log10 ratio of the January and February Kv/Av (b over a) in depth
bands, from identical starts (crunveg 2119-12-31, old build).  Known: 1800 s vs 1200 s
gives -0.015 to -0.021 below 100 m in three independent pairs; noise +-0.002.
A change that removes the step dependence brings its 1800/1200 pair to ~0.
"""
import os, sys
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4'
PAIRS = [('base: 1200 -> 1800 s (iter 5)', 'TEST_dt1200', 'TEST_dt1800'),
         ('1200 iter5 -> 1800 iter8', 'TEST_dt1200', 'TEST_dt1800h8'),
         ('horiz. off: 1200 -> 1800 s (iter 0)', 'TEST_dt1200h0', 'TEST_dt1800h0'),
         ('effect of horiz. off at 1200 s', 'TEST_dt1200', 'TEST_dt1200h0'),
         ('effect of horiz. off at 1800 s', 'TEST_dt1800', 'TEST_dt1800h0')]
# Patched build (FESOM fix/idemix-hor-smoothing-dt 8bb5f16a on 766020ec): `step_screens.py xd`
PAIRS_XD = [('patched build, legacy settings vs old build (1200 s): must be ~0', 'TEST_dt1200', 'TEST_xd1200'),
            ('FIX: 1200 legacy -> 1800 smooth_dtref 1200', 'TEST_dt1200', 'TEST_xd1800s'),
            ('FIX + halo: 1200 halo -> 1800 smooth_dtref 1200 + halo', 'TEST_xd1200x', 'TEST_xd1800sx'),
            ('effect of halo exchange at 1200 s', 'TEST_dt1200', 'TEST_xd1200x'),
            ('effect of removing the smoothing at 1800 s', 'TEST_dt1800', 'TEST_xd1800s0'),
            ('base again: 1200 -> 1800 s legacy', 'TEST_dt1200', 'TEST_dt1800')]
if sys.argv[1:] == ['xd']:
    PAIRS = PAIRS_XD
BANDS = ((0, 100), (100, 500), (500, 1500), (1500, 7000))


def load(exp, v):
    with xr.open_dataset(f'{R}/{exp}/outdata/fesom/{v}.fesom.2120.nc', decode_times=False) as d:
        a = d[v].values.astype('f8'); dz = [c for c in d[v].dims if c.startswith('nz')][0]
        dep = np.abs(d[dz].values)
    if a.shape[1] == len(dep): a = np.swapaxes(a, 1, 2)
    return np.where(np.abs(a) > 1e30, np.nan, a), dep


for lab, ea, eb in PAIRS:
    print(lab)
    if 'must be ~0' in lab:
        a, _ = load(ea, 'Kv'); b, _ = load(eb, 'Kv')
        print(f'   Kv bit-identical: {np.array_equal(np.nan_to_num(a), np.nan_to_num(b))}')
    for v in ('Kv', 'Av'):
        a, dep = load(ea, v); b, _ = load(eb, v)
        for m, mn in ((0, 'Jan'), (1, 'Feb')):
            out = []
            for z0, z1 in BANDS:
                kz = (dep >= z0) & (dep < z1)
                x, y = a[m][:, kz], b[m][:, kz]
                ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
                out.append(f'{z0}-{z1}m {np.mean(np.log10(y[ok] / x[ok])):+7.4f}')
            print(f'   {v} {mn}  ' + '  '.join(out), flush=True)
