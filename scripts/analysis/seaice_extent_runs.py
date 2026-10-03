"""Sea-ice extent guards for a set of runs: NH March and SH September, per year.

Extent = area of cells with a_ice >= 0.15 from FESOM's regular daily output (cos-lat cell
areas, ocean cells only), monthly means of the daily field.  Prints per-year values and the
window means, branch minus control.  Protocol references (OSI-SAF 1990-2019): NH March
about 15.0 M km2 on this measure was the campaign's working number; compare within the
table first.

Usage:  CONTROL=PI200 RUNS=PI200_gmR2500,PI200_gmR1500 Y0=1580 Y1=1599 python3 scripts/analysis/seaice_extent_runs.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
CONTROL = os.environ.get('CONTROL', 'PI200'); RUNS = os.environ.get('RUNS', 'PI200_gmR2500,PI200_gmR1500').split(',')
Y0, Y1 = int(os.environ.get('Y0', 1580)), int(os.environ.get('Y1', 1599))


def extents(arm, y):
    with xr.open_dataset(f'{R}/{arm}/outdata/fesom/a_ice.fesom.gr.{y}.nc', decode_times=False) as d:
        a = d['a_ice'].values.astype('f8'); lat = d['lat'].values
    a[np.abs(a) > 100] = np.nan
    n = [31, 29 if a.shape[0] == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; e = np.cumsum([0] + n)
    mar = np.nanmean(a[e[2]:e[3]], 0); sep = np.nanmean(a[e[8]:e[9]], 0)
    A = np.cos(np.deg2rad(lat))[:, None] * (0.5 * 111.195e3) ** 2 * np.ones(a.shape[2])[None, :]
    nh = lat[:, None] > 0
    return float(np.nansum(A[(mar >= 0.15) & nh]) / 1e12), float(np.nansum(A[(sep >= 0.15) & ~nh]) / 1e12)


D = {a: {y: extents(a, y) for y in range(Y0, Y1 + 1) if os.path.exists(f'{R}/{a}/outdata/fesom/a_ice.fesom.gr.{y}.nc')} for a in [CONTROL] + RUNS}
print(f'{Y0}-{Y1}: NH March / SH September extent [M km2]')
print(f'  {"year":<6}' + ''.join(f'{a:>22}' for a in D))
for y in range(Y0, Y1 + 1):
    print(f'  {y:<6}' + ''.join(f'{D[a][y][0]:11.2f}{D[a][y][1]:11.2f}' if y in D[a] else f'{"":>22}' for a in D))
for lab, yy in ((f'{Y0}-{Y1}', range(Y0, Y1 + 1)), (f'{Y0 + 10}-{Y1}', range(Y0 + 10, Y1 + 1))):
    c = np.array([D[CONTROL][y] for y in yy if y in D[CONTROL]])
    print(f'  mean {lab}: {CONTROL} NH Mar {c[:, 0].mean():.2f} SH Sep {c[:, 1].mean():.2f}; ' + '; '.join(
        f'{a} {np.array([D[a][y] for y in yy if y in D[a]])[:, 0].mean() - c[:, 0].mean():+.2f} / {np.array([D[a][y] for y in yy if y in D[a]])[:, 1].mean() - c[:, 1].mean():+.2f}' for a in RUNS)
          + f'   (ctl sd {c[:, 0].std(ddof=1):.2f} / {c[:, 1].std(ddof=1):.2f})')
