"""Ocean temperature drift by depth band along the production line, grouped by vertical mixing scheme.

The production line ran four mixing configurations in sequence, which makes a long record for each:
  KPP (FESOM's own)        1850-2059   (PICAL, PICAL_momixoff, PICAL_ccnice)
  cvmix_TKE                2060-2079   (PICAL_ccnice)
  cvmix_TKE + IDEMIX       2080-2139   (PICAL_ccnice, PICAL_crunveg, PICAL_crunveg_gmhemi1800)
  cvmix_TKE                2140-2209   (PICAL_crunveg_gmhemi_tke, PICAL_crunveg_tke_albsn082)
  cvmix_TKE + IDEMIX       2170-2209   (PICAL_crunveg_idemix_albsn082, branch of the line above)
Other things changed inside these eras as well (sea-ice CCN 1940, snow albedo 1970/1980, h0min
2000-2019, aerosol 2020, albedo/land 2100, hemispheric GM and 1800 s step 2120, albsn 0.82 2160),
so an era trend is not a clean attribution to the mixing scheme; the paired 20-year arms
(mixing_arms_temp_trend_by_depth.txt) are the clean test, this is the long view.

Global-mean temperature profile per year from reval's spin-up cache, averaged over each band by layer
thickness. Trend in mK/decade +- 2 s.e. (no autocorrelation correction); and the band mean at the end of
each segment minus the first decade of the record (1850-1859), i.e. the accumulated drift, in mK.

Usage: mixing_eras_temp_drift_by_depth.py   (reval environment; writes data/clim/mixing_eras_temp_drift_by_depth.txt)
"""
import os
import numpy as np
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
C = f'{P}/reval/spinup_cache'
BANDS = [(0, 100), (100, 300), (300, 700), (700, 1500), (1500, 2500), (2500, 6000)]
ZI = np.abs(xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')['nz'].values)
DZ, ZM = np.diff(ZI), 0.5 * (ZI[1:] + ZI[:-1])
SEG = [('KPP 1860-1939', 'critical_path_tke', 1860, 1939),
       ('KPP 1940-1999', 'critical_path_tke', 1940, 1999),
       ('KPP 2000-2059', 'critical_path_tke', 2000, 2059),
       ('KPP 2020-2059 (settled)', 'critical_path_tke', 2020, 2059),
       ('TKE 2060-2079', 'critical_path_tke', 2060, 2079),
       ('TKE+IDEMIX 2080-2099', 'critical_path_tke', 2080, 2099),
       ('TKE+IDEMIX 2100-2139', 'critical_path_tke', 2100, 2139),
       ('TKE 2140-2159', 'critical_path_tke', 2140, 2159),
       ('TKE 2160-2209', 'critical_path_tke', 2160, 2209),
       ('TKE+IDEMIX 2170-2209 (arm)', 'critical_path_idemix', 2170, 2209)]


def bands(cache, y):
    prof = np.load(f'{C}/{cache}/hovm_temp_{y}.npy')
    n = len(prof)
    return [float(np.sum(prof[(ZM[:n] >= a) & (ZM[:n] < b)] * DZ[:n][(ZM[:n] >= a) & (ZM[:n] < b)]) /
                  np.sum(DZ[:n][(ZM[:n] >= a) & (ZM[:n] < b)])) for a, b in BANDS]


def trend(y, x):
    A = np.vstack([y - y.mean(), np.ones_like(y)]).T
    c, res = np.linalg.lstsq(A, x, rcond=None)[:2]
    return 10 * c[0], 20 * np.sqrt(res[0] / (len(y) - 2) / np.sum((y - y.mean()) ** 2))


ref = np.mean([bands('critical_path_tke', y) for y in range(1850, 1860)], axis=0)
hdr = f"{'segment':28s}" + ''.join(f"{f'{a}-{b} m':>15s}" for a, b in BANDS)
out = ['Global-mean ocean temperature by depth band along the production line, by mixing era.', '',
       'Trend, mK/decade +- 2 s.e.', hdr]
ends = []
for name, cache, y0, y1 in SEG:
    yy = np.arange(y0, y1 + 1, dtype=float)
    T = np.array([bands(cache, int(y)) for y in yy])
    out.append(f'{name:28s}' + ''.join('{:>+9.1f}±{:<5.1f}'.format(*[1e3 * v for v in trend(yy, T[:, k])]) for k in range(len(BANDS))))
    ends.append((name, 1e3 * (T[-10:].mean(0) - ref)))
out += ['', 'Last 10 years of the segment minus 1850-1859, mK (accumulated drift since the start of the record)', hdr]
for name, d in ends:
    out.append(f'{name:28s}' + ''.join(f'{v:>+15.0f}' for v in d))
out += ['', '1850-1859 band means, C', f"{'':28s}" + ''.join(f'{v:>15.3f}' for v in ref)]
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/mixing_eras_temp_drift_by_depth.txt', 'w').write(txt + '\n')
