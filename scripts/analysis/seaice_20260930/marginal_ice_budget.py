"""Is the extra marginal-sea winter ice after the 2120 branch imported or frozen locally?

FESOM's own ice-volume tendencies split the answer: dyngrice (advection/divergence, i.e.
ice brought in by the drift) and thdgrice (thermodynamic growth/melt in place), both
grid-mean volume per unit area.  Summed over Nov-Apr for boxes where ice expanded
(Labrador interior, Labrador+Newfoundland, Bering south of 66N) and a central-Arctic
reference, per run.  Also the mean ice speed (daily |u_ice|) and concentration.

Usage: python marginal_ice_budget.py            (runs and years set below)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4'
RUNS = [('PICAL_crunveg', 2110, 2119, 'old build 1200 s'),
        ('PICAL_crunveg_ihf0', 2121, 2129, 'new build 1800 s'),
        ('PICAL_crunveg_ihf1', 2121, 2129, 'new build 1800 s + McPhee')]
BOX = {'Labrador interior': (56, 62, 300, 310), 'Labrador+Newfoundland': (40, 66, 290, 316),
       'Bering <66N': (50, 66, 163, 203), 'central Arctic >80N': (80, 90, 0, 360)}
WIN = [10, 11, 0, 1, 2, 3]          # Nov-Apr
SEC = 86400 * 30.44


def rd(exp, v, y):
    with xr.open_dataset(f'{R}/{exp}/outdata/fesom/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
        a = np.squeeze(d[v].values).astype('f8'); lat = d['lat'].values; lon = d['lon'].values % 360
    a = np.where(np.abs(a) > 1e30, np.nan, a)
    if a.shape[0] != 12:
        ml = [31, 29 if a.shape[0] == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        e = np.cumsum([0] + ml); a = np.stack([np.nanmean(a[e[i]:e[i + 1]], 0) for i in range(12)])
    return a, lat, lon


for exp, y0, y1, lab in RUNS:
    acc = {}
    for y in range(y0, y1 + 1):
        for v in ('thdgrice', 'dyngrice', 'a_ice'):
            a, lat, lon = rd(exp, v, y); acc[v] = acc.get(v, 0) + a / (y1 - y0 + 1)
        u, _, _ = rd(exp, 'uice', y); w, _, _ = rd(exp, 'vice', y)
        acc['spd'] = acc.get('spd', 0) + np.hypot(u, w) / (y1 - y0 + 1)   # monthly means of daily components
    LAT, LON = np.meshgrid(lat, lon, indexing='ij'); W = np.cos(np.deg2rad(LAT))
    print(f'{exp} {y0}-{y1} ({lab})')
    for name, (la0, la1, lo0, lo1) in BOX.items():
        k = (LAT >= la0) & (LAT <= la1) & (LON >= lo0) & (LON <= lo1)
        f = lambda v, m: np.nansum((np.nan_to_num(acc[v][m]) * W)[k]) / W[k].sum()
        thd = sum(f('thdgrice', m) for m in WIN) * SEC * 100
        dyn = sum(f('dyngrice', m) for m in WIN) * SEC * 100
        ai = np.mean([f('a_ice', m) for m in (1, 2)])
        kk = k & (acc['a_ice'][2] > 0.15)
        spd = np.nansum((acc['spd'][2] * W)[kk]) / max(W[kk].sum(), 1e-9) * 100
        print(f'   {name:22s} Nov-Apr thermo {thd:7.1f} cm  dyn {dyn:7.1f} cm   Feb-Mar conc {ai:5.2f}   Mar ice speed {spd:5.1f} cm/s')
    print()
