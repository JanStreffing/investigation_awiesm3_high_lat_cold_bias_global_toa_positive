"""Ice drift in the 2-month albedo tests: does the time step or the FESOM build slow the pack?

All tests start from the same 2100-01-01 restart; ice drift responds within days, so days
15-59 of Jan-Feb 2100 separate the two factors: old build (libfesom 4c76f03d) at 1200 s
(perfE, perfG) and 1800 s (perfF, perfH), new build (18087062) at 1800 s (nolockB1/B4,
best25).  Mean of daily |u_ice| over a_ice >= 0.8, per region.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4'
T = [('TEST_perfE', 'old 1200'), ('TEST_perfG', 'old 1200'), ('TEST_perfF', 'old 1800'),
     ('TEST_perfH', 'old 1800'), ('TEST_nolockB1', 'new 1800'), ('TEST_nolockB4', 'new 1800'),
     ('TEST_best25', 'new 1800')]
BOX = {'>80N': (80, 90, 0, 360), 'Beaufort': (70, 80, 190, 235), 'Labrador': (52, 66, 290, 316),
       'Bering': (55, 66, 163, 203), 'Weddell': (-78, -60, 300, 360)}
for t, lab in T:
    def rd(v):
        with xr.open_dataset(f'{R}/{t}/outdata/fesom/{v}.fesom.gr.2100.nc', decode_times=False) as d:
            a = np.squeeze(d[v].values).astype('f8'); return np.where(np.abs(a) > 1e30, np.nan, a), d['lat'].values, d['lon'].values % 360
    u, lat, lon = rd('uice'); v, _, _ = rd('vice'); a, _, _ = rd('a_ice')
    n = min(u.shape[0], a.shape[0], 59)
    s = np.hypot(u[14:n], v[14:n]); aa = a[14:n]
    LAT, LON = np.meshgrid(lat, lon, indexing='ij'); W = np.cos(np.deg2rad(LAT))
    out = []
    for name, (la0, la1, lo0, lo1) in BOX.items():
        k = (LAT >= la0) & (LAT <= la1) & (LON >= lo0) & (LON <= lo1)
        m = (aa >= 0.8) & k[None]
        out.append(f'{name} {100 * np.nansum(np.where(m, s, 0) * W[None]) / np.sum(m * W[None]):5.2f}')
    print(f'{t:14s} {lab:9s} days 15-{n}  |u_ice| cm/s: ' + '  '.join(out), flush=True)
