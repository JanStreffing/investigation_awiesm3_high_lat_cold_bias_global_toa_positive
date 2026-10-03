"""Monthly Labrador sea ice (and MLD2) for any FESOM mesh, for the run-history search.

Regions from node lon/lat: labsea 52-66N 65-45W, interior 56-62N 60-50W.
Weights: polygon area of each node's dual cell from bounds_lon/bounds_lat (local
equirectangular, fine at box scale); uniform if the file has no bounds. Coordinates
come from the data file itself or, for files written without them, from a donor file
on the same mesh (same nod2). a_ice: masked = no ice. MLD2: masked excluded, |MLD|.
Usage: python labrador_any_mesh.py <coord_file|-> <out.csv> <exp> <dir> <y0> <y1> [stride]
"""
import sys, calendar, numpy as np
from netCDF4 import Dataset
cf, outf, exp, d, y0, y1 = sys.argv[1:7]; y0, y1 = int(y0), int(y1)
stride = int(sys.argv[7]) if len(sys.argv) > 7 else 1
def coords(path):
    with Dataset(path) as nc:
        lon = nc.variables['lon'][:].astype('f8'); lat = nc.variables['lat'][:].astype('f8')
        bl = nc.variables['bounds_lon'][:].astype('f8') if 'bounds_lon' in nc.variables else None
        bt = nc.variables['bounds_lat'][:].astype('f8') if 'bounds_lat' in nc.variables else None
    return np.asarray(lon), np.asarray(lat), bl, bt
src = cf if cf != '-' else f'{d}/a_ice.fesom.{y0}.nc'
lon, lat, bl, bt = coords(src); lon = ((lon + 180) % 360) - 180
R = {'labsea': (lat >= 52) & (lat <= 66) & (lon >= -65) & (lon <= -45),
     'interior': (lat >= 56) & (lat < 62) & (lon >= -60) & (lon < -50)}
W = {}
for r, m in R.items():
    idx = np.where(m)[0]
    if bl is None:
        w = np.ones(len(idx))
    else:
        x = np.asarray(bl)[idx]; y = np.asarray(bt)[idx]
        bad = ~np.isfinite(x) | (np.abs(x) > 1e5) | ~np.isfinite(y) | (np.abs(y) > 1e5)
        for k in range(1, x.shape[1]):   # fill padding with the previous vertex
            x[:, k] = np.where(bad[:, k], x[:, k-1], x[:, k]); y[:, k] = np.where(bad[:, k], y[:, k-1], y[:, k])
        x = (((x - lon[idx, None]) + 180) % 360) - 180   # relative lon, no dateline jumps
        xc = np.deg2rad(x) * np.cos(np.deg2rad(lat[idx, None])); yc = np.deg2rad(y)
        w = 0.5 * np.abs(np.sum(xc * np.roll(yc, -1, 1) - np.roll(xc, -1, 1) * yc, axis=1)) * 6.371e6**2
    W[r] = (idx, w)
out = open(outf, 'a')
for y in range(y0, y1 + 1, stride):
    for var in ('a_ice', 'MLD2'):
        try:
            with Dataset(f'{d}/{var}.fesom.{y}.nc') as nc:
                v = nc.variables[var][:]
                a = v.filled(np.nan) if hasattr(v, 'filled') else np.asarray(v)
        except (OSError, KeyError):
            continue
        a = np.asarray(a, dtype='f8')
        if a.ndim == 2 and a.shape[0] != len(lon):   # (time, nod2)
            pass
        nt = a.shape[0]
        if nt in (365, 366):
            ml = [calendar.monthrange(2001 if nt == 365 else 2000, m)[1] for m in range(1, 13)]
            e = np.cumsum([0] + ml)
        elif nt == 12:
            e = np.arange(13)
        else:
            print(f'{exp} {y} {var} nt={nt} skipped', flush=True); continue
        for r, (idx, w) in W.items():
            x = a[:, idx]
            if var == 'a_ice':
                x = np.nan_to_num(x, nan=0.0); x[np.abs(x) > 1e10] = 0.0
                s = x.dot(w) / w.sum()
            else:
                x[np.abs(x) > 1e10] = np.nan; x = np.abs(x); ok = np.isfinite(x)
                s = np.where(ok, x, 0).dot(w) / np.maximum((ok * w).sum(axis=1), 1e-30)
            for m in range(12):
                out.write(f'{exp},{y},{m+1},{r},{var},{np.mean(s[e[m]:e[m+1]]):.5f}\n')
    out.flush(); print(f'{exp} {y} ok', flush=True)
