"""Southern Ocean winter state from FESOM node output, for runs of any mesh: is the pack venting?

Over Y0-Y1, 60-78S open-ocean nodes (no cavity), area-weighted:
  JJA mixed layer depth (MLD2) under the pack (JJA a_ice >= 0.8) and under all ice (>= 0.15);
  JJA floe thickness m_ice/a_ice and snow m_snow/a_ice under ice; September ice area and extent;
  annual 0-50 m salinity and temperature and 100-500 m temperature against PHC3 on the same mesh.
Node files: a_ice, m_ice, m_snow (daily or monthly), MLD2 (monthly), temp, salt (monthly).

Usage:  ROOT=... ARM=PI Y0=2000 Y1=2014 MESH=<mesh diag> PHC_DIR=<climatology dir> python3 scripts/analysis/so_winter_state_nodes.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ['ROOT']; ARM = os.environ['ARM']; MESH = os.environ['MESH']; PHC = os.environ['PHC_DIR']
Y0, Y1 = int(os.environ['Y0']), int(os.environ['Y1'])
with xr.open_dataset(MESH) as m:
    lat = m['lat'].values; area = m['nod_area'].values[0]; nlv = m['nlevels_nod2D'].values; z = m['nz1'].values
    ul = m['ulevels_nod2D'].values if 'ulevels_nod2D' in m else np.ones(len(lat), int)
NZ = len(z); k = (lat >= -78) & (lat <= -60) & (ul == 1); N = Y1 - Y0 + 1
JJA, SEP = [5, 6, 7], [8]


def rd(v, y, months):
    with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/{v}.fesom.{y}.nc', decode_times=False) as d:
        a = d[v]
        dims = a.dims; a = a.values.astype('f8')
    a[np.abs(a) > 1e10] = np.nan
    if a.shape[0] > 12:                                           # daily
        n = [31, 29 if a.shape[0] == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; e = np.cumsum([0] + n)
        a = np.stack([np.nanmean(a[e[i]:e[i + 1]], 0) for i in range(12)])
    a = np.nanmean(a[months], 0)
    if a.ndim == 2 and a.shape[0] != len(lat): a = a.T                # (nod2, nz)
    return a


def col(v, y):
    a = rd(v, y, list(range(12)))
    return np.where(np.arange(a.shape[1])[None, :] < (nlv - 1)[:, None], a[:, :NZ], np.nan)


def lay(a, lo, hi):
    kz = (z >= lo) & (z < hi); dz = np.diff(np.concatenate([[0.0], 0.5 * (z[1:] + z[:-1]), [z[-1] + 0.5 * (z[-1] - z[-2])]]))
    w = area[k][:, None] * dz[None, kz] * np.isfinite(a[k][:, kz]); return np.nansum(np.nan_to_num(a[k][:, kz]) * w) / w.sum()


A = {'a': 0, 'h': 0, 's': 0, 'mld': 0, 'asep': 0, 'T': 0, 'S': 0}
for y in range(Y0, Y1 + 1):
    A['a'] = A['a'] + rd('a_ice', y, JJA) / N; A['h'] = A['h'] + rd('m_ice', y, JJA) / N; A['s'] = A['s'] + rd('m_snow', y, JJA) / N
    A['mld'] = A['mld'] + np.abs(rd('MLD2', y, JJA)) / N; A['asep'] = A['asep'] + rd('a_ice', y, SEP) / N
    A['T'] = A['T'] + col('temp', y) / N; A['S'] = A['S'] + col('salt', y) / N
P = {}
for v, n in (('temp', 'T'), ('salt', 'S')):
    with xr.open_dataset(f'{PHC}/{v}.fesom.1958.nc', decode_times=False) as d:
        p = np.squeeze(d[v].values).astype('f8')
    if p.shape[0] != len(lat): p = p.T
    p = p[:, :NZ]; p[np.abs(p) > 100] = np.nan; P[n] = np.where(np.arange(NZ)[None, :] < (nlv - 1)[:, None], p, np.nan)
av = lambda f, kk: float(np.average(f[kk], weights=area[kk]))
pack = k & (A['a'] >= 0.8); ice = k & (A['a'] >= 0.15)
print(f'{ARM} {Y0}-{Y1}, 60-78S open ocean, {len(lat)} nodes, spacing {np.sqrt(np.average(area[k], weights=area[k])) / 1e3:.0f} km')
print(f'  JJA: pack area {area[pack].sum() / 1e12:.2f} M km2, MLD under pack {av(A["mld"], pack):.0f} m, under all ice {av(A["mld"], ice):.0f} m, '
      f'floe thickness pack {av(A["h"] / np.maximum(A["a"], 0.15), pack):.2f} m, snow {av(A["s"] / np.maximum(A["a"], 0.15), pack):.2f} m, mean conc. under ice {av(A["a"], ice):.2f}')
print(f'  Sep: SH ice area {np.nansum(A["asep"][k] * area[k]) / 1e12:.2f}, extent {area[k & (A["asep"] >= 0.15)].sum() / 1e12:.2f} M km2 (60-78S only)')
print(f'  column vs PHC3: T 0-50 {lay(A["T"], 0, 50):+.2f} ({lay(A["T"], 0, 50) - lay(P["T"], 0, 50):+.2f}), T 100-500 {lay(A["T"], 100, 500):+.2f} ({lay(A["T"], 100, 500) - lay(P["T"], 100, 500):+.2f}), '
      f'T 500-1500 {lay(A["T"], 500, 1500):+.2f} ({lay(A["T"], 500, 1500) - lay(P["T"], 500, 1500):+.2f}), S 0-50 {lay(A["S"], 0, 50):.3f} ({lay(A["S"], 0, 50) - lay(P["S"], 0, 50):+.3f}), '
      f'S 200-500 {lay(A["S"], 200, 500):.3f} ({lay(A["S"], 200, 500) - lay(P["S"], 200, 500):+.3f})')
