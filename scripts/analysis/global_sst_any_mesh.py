"""Annual global-mean SST for any FESOM mesh: the global-temperature control for the
Labrador sea-ice history search (labrador_any_mesh.py). Area weights from the dual-cell
bounds polygons, coordinates from a donor file when the data file has none. Masked
(cavity) values are excluded. Also the 50-90N and 50-90S means.
Usage: python global_sst_any_mesh.py <coord_file|mesh.nc|-> <out.csv> <exp> <dir> <y0> <y1> [stride]
"""
import sys, numpy as np
from netCDF4 import Dataset
cf, outf, exp, d, y0, y1 = sys.argv[1:7]; y0, y1 = int(y0), int(y1)
stride = int(sys.argv[7]) if len(sys.argv) > 7 else 1
src = cf if cf != '-' else f'{d}/sst.fesom.{y0}.nc'
with Dataset(src) as nc:
    lon = np.asarray(nc.variables['lon'][:], 'f8'); lat = np.asarray(nc.variables['lat'][:], 'f8')
    area = np.asarray(nc.variables['cell_area'][:], 'f8') if 'cell_area' in nc.variables else None
    if area is None:
        bl = np.asarray(nc.variables['bounds_lon'][:], 'f8'); bt = np.asarray(nc.variables['bounds_lat'][:], 'f8')
if area is not None:   # mesh.nc (CDO grid description) carries the dual-cell areas
    lon = np.rad2deg(lon) if np.abs(lon).max() < 7 else lon; lat = np.rad2deg(lat) if np.abs(lat).max() < 1.6 else lat
    bl = bt = None
if bl is not None:
  bad = ~np.isfinite(bl) | (np.abs(bl) > 1e5) | ~np.isfinite(bt) | (np.abs(bt) > 1e5)
  for k in range(1, bl.shape[1]):
      bl[:, k] = np.where(bad[:, k], bl[:, k-1], bl[:, k]); bt[:, k] = np.where(bad[:, k], bt[:, k-1], bt[:, k])
  x = ((bl - lon[:, None] + 180) % 360) - 180
  xc = np.deg2rad(x) * np.cos(np.deg2rad(lat[:, None])); yc = np.deg2rad(bt)
  w = 0.5 * np.abs(np.sum(xc * np.roll(yc, -1, 1) - np.roll(xc, -1, 1) * yc, axis=1))
  w[~np.isfinite(w) | (np.abs(lat) > 89.5)] = 0.0   # polar-cap polygons are unreliable; tiny area
else:
    w = area.copy()
M = {'global': np.ones_like(lat, bool), 'nh50': lat >= 50, 'sh50': lat <= -50}
out = open(outf, 'a')
for y in range(y0, y1 + 1, stride):
    try:
        with Dataset(f'{d}/sst.fesom.{y}.nc') as nc:
            v = nc.variables['sst'][:]
            a = np.asarray(v.filled(np.nan) if hasattr(v, 'filled') else v, 'f8')
    except (OSError, KeyError):
        continue
    a[np.abs(a) > 1e10] = np.nan
    am = np.nanmean(a, axis=0)
    ok = np.isfinite(am)
    vals = [np.sum(np.where(ok & m, am, 0) * w) / np.sum(w * (ok & m)) for m in M.values()]
    out.write(f'{exp},{y},' + ','.join(f'{v:.4f}' for v in vals) + '\n'); out.flush()
    print(f'{exp} {y} ok', flush=True)
