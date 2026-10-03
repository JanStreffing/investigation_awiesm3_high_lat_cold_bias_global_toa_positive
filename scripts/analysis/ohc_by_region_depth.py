"""Annual ocean heat content by region and depth band, plus regional net surface heat flux,
to locate where GM 2500 reduces ocean heat uptake (CORE3 mesh, 220509 nodes).

Regions (node lon/lat): SO < 40S; SPNA 45-66N 80W-10E (Labrador, Irminger, Iceland basin);
NORDARC > 66N (Nordic Seas + Arctic); REST everything else. Depth bands 0-700, 700-2000,
> 2000 m by layer mid-depth. Volume = hnode (time mean of the first year) x nod_area of the
layer's upper interface. OHC = 4.1e6 J m-3 K-1 x T x V. fh: net surface heat flux into the
ocean, area-weighted with the surface nod_area.
Usage: python ohc_by_region_depth.py <mesh.diag.nc> <out.csv> <exp> <outdata_fesom> <y0> <y1>
Output rows: exp,year,kind,region,band,value   (kind ohc [J] or fh [W])
"""
import sys, numpy as np
from netCDF4 import Dataset
diag, outf, exp, d, y0, y1 = sys.argv[1:7]; y0, y1 = int(y0), int(y1)
with Dataset(diag) as nc:
    area = np.asarray(nc.variables['nod_area'][:], 'f8')          # (nz=48, nod2)
    lon = np.asarray(nc.variables['lon'][:], 'f8'); lat = np.asarray(nc.variables['lat'][:], 'f8')
    zb = np.abs(np.asarray(nc.variables['nz'][:], 'f8'))          # interface depths (48)
if np.abs(lon).max() < 7: lon, lat = np.rad2deg(lon), np.rad2deg(lat)
lon = ((lon + 180) % 360) - 180
with Dataset(f'{d}/hnode.fesom.{y0}.nc') as nc:
    h = np.asarray(nc.variables['hnode'][:], 'f8').mean(0)       # (nz1, nod2) or (nod2, nz1)
if h.shape[0] != 47: h = h.T
h = np.where(np.isfinite(h) & (np.abs(h) < 1e5), h, 0.0)
vol = h * area[:47]
zmid = 0.5 * (zb[:-1] + zb[1:])
REG = {'SO': lat < -40, 'SPNA': (lat >= 45) & (lat < 66) & (lon >= -80) & (lon <= 10), 'NORDARC': lat >= 66}
REG['REST'] = ~(REG['SO'] | REG['SPNA'] | REG['NORDARC'])
BAND = {'0-700': zmid < 700, '700-2000': (zmid >= 700) & (zmid < 2000), '2000+': zmid >= 2000}
out = open(outf, 'a')
for y in range(y0, y1 + 1):
    try:
        with Dataset(f'{d}/temp.fesom.{y}.nc') as nc:
            v = nc.variables['temp']
            t = np.zeros(v.shape[1:], 'f8')
            for m in range(v.shape[0]):
                x = np.asarray(v[m], 'f8'); x[~np.isfinite(x) | (np.abs(x) > 1e5)] = 0.0; t += x
            t /= v.shape[0]
        if t.shape[0] != 47: t = t.T
        with Dataset(f'{d}/fh.fesom.{y}.nc') as nc:
            f = np.asarray(nc.variables['fh'][:], 'f8'); f[~np.isfinite(f) | (np.abs(f) > 1e5)] = 0.0; f = f.mean(0)
    except (OSError, KeyError) as e:
        print(f'{exp} {y} {e}', flush=True); continue
    hc = 4.1e6 * t * vol
    for r, m in REG.items():
        for b, k in BAND.items():
            out.write(f'{exp},{y},ohc,{r},{b},{hc[k][:, m].sum():.6e}\n')
        out.write(f'{exp},{y},fh,{r},sfc,{(f[m] * area[0][m]).sum():.6e}\n')
    out.flush(); print(f'{exp} {y} ok', flush=True)
