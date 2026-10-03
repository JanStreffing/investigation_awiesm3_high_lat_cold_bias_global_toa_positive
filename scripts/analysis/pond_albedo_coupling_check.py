"""Does OpenIFS see FESOM's melt ponds?  July pond fraction (FESOM) against surface albedo (OpenIFS).

In the coupled build FESOM's ice growth uses the atmosphere flux directly
(ice_thermo_cpl.F90: Qatmice = -a2ihf), so FESOM's pond albedo can only matter if OpenIFS
uses the A_Ice_albedo FESOM sends (radpar.F90: IF (LNEMOLIMALB)).  This compares, per
decade, over Arctic pack ice in July:
  apnd  FESOM melt-pond area fraction, nodes with a_ice > 0.8 and lat > 70N, cell_area weighted
  fal   OpenIFS surface albedo, cells with ci > 0.8 and lat > 70N, cos-lat weighted
If apnd moves when the pond settings change but fal does not, the atmosphere is not using
FESOM's albedo.  16E_1990 changed pond settings at 1400.

Usage:  python3 scripts/analysis/pond_albedo_coupling_check.py
        RUN=16E_1990 DECADES=1390,1400,1410 python3 ...
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
RUN = os.environ.get('RUN', '16E_1990')
DECADES = [int(d) for d in os.environ.get('DECADES', '1390,1400,1410').split(',')]
with xr.open_dataset('/work/ab0246/a270092/input/fesom2/core3/mesh.nc') as m:
    AREA = m['cell_area'].values.astype('f8')
JUL = 6


def fesom_july(var, y):
    with xr.open_dataset(f'{R}/{RUN}/outdata/fesom/{var}.fesom.{y}.nc', decode_times=False) as d:
        v = np.squeeze(d[var].values)
        lat = d['lat'].values
    n = v.shape[0]
    if n == 12:                       # monthly file (apnd)
        return v[JUL], lat
    j0 = 181 + (1 if n == 366 else 0)  # daily file (a_ice): July = days 182-212 (183-213 in a leap year)
    return v[j0:j0 + 31].mean(0), lat


def oifs_july(var, y):
    with xr.open_dataset(f'{R}/{RUN}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc', decode_times=False) as d:
        k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values)[JUL], np.squeeze(d['lat'].values)


print(f'{RUN}: July, Arctic pack ice (>0.8 concentration, north of 70N)')
print(f'{"decade":>10s} {"FESOM apnd":>11s} {"OIFS fal":>9s}')
for d0 in DECADES:
    ap, fa = [], []
    for y in range(d0, d0 + 10):
        a, lat = fesom_july('a_ice', y); p, _ = fesom_july('apnd', y)
        k = (a > 0.8) & (lat > 70) & np.isfinite(p)
        w = AREA[:a.size] if AREA.size == a.size else np.ones(a.size)
        ap.append(np.sum(p[k] * w[k]) / np.sum(w[k]))
        ci, la = oifs_july('ci', y); fal, _ = oifs_july('fal', y)
        L = np.broadcast_to(la[:, None], ci.shape)
        kk = (ci > 0.8) & (L > 70) & np.isfinite(fal)
        fa.append(np.average(fal[kk], weights=np.cos(np.deg2rad(L[kk]))))
    print(f'{d0}-{d0 + 9} {np.mean(ap):11.4f} {np.mean(fa):9.4f}   (year-to-year sd {np.std(ap):.4f} / {np.std(fa):.4f})')
