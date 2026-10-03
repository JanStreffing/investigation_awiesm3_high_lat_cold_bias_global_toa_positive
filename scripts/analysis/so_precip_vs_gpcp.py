"""Model precipitation over the Southern Ocean against GPCP, by band, ocean only.

Model lsp+cp on the remapped grid (accumulated m per hourly step unless ACC=1), annual mean
over Y0-Y1; GPCP v2.3 monthly climatology 1990-2014 interpolated to the model grid.  Prints
mm/yr and Gt/yr over the ocean in 60-78S, 45-60S and south of 60S including Antarctic land
(where GPCP is least reliable).

Usage:  ROOT=... ARM=PI200 Y0=1585 Y1=1599 python3 scripts/analysis/so_precip_vs_gpcp.py
"""
import os, glob
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ['ROOT']; ARM = os.environ['ARM']; Y0, Y1 = int(os.environ['Y0']), int(os.environ['Y1']); ACC = float(os.environ.get('ACC', 3600.0))
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'


def ann(v):
    acc = 0; n = 0
    for y in range(Y0, Y1 + 1):
        p = glob.glob(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc') + glob.glob(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc')
        if not p: continue
        with xr.open_dataset(p[0], decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8'); lat = np.squeeze(d['lat'].values); lon = np.squeeze(d['lon'].values)
        acc = acc + a.mean(0) / ACC; n += 1
    return acc / n, lat, lon


P, lat, lon = ann('lsp'); C, _, _ = ann('cp'); P = (P + C) * 365.25 * 86400 * 1000            # mm/yr
o = np.argsort(lon % 360); P = P[:, o]; lon = (lon % 360)[o]
with xr.open_dataset('/work/ab0246/a270092/obs/gpcp/precip.mon.mean.nc') as d:
    g = d['precip'].sel(time=d['time.year'].isin(range(1990, 2015))).mean('time')
    g = g.assign_coords(lon=g['lon'] % 360).sortby('lon').sortby('lat')
    G = np.asarray(g.interp(lat=lat, lon=lon).values, float) * 365.25                       # mm/day -> mm/yr
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = (lsm[0] if lsm.ndim == 3 else lsm)[:, o]
dlat = np.abs(np.gradient(lat)); A = (np.cos(np.deg2rad(lat)) * dlat * 111.195e3 * (360.0 / len(lon)) * 111.195e3)[:, None] * np.ones(len(lon))[None, :]
print(f'{ARM} {Y0}-{Y1}: precipitation against GPCP 1990-2014')
for name, k in (('ocean 60-78S', (lsm <= 0.5) & (lat[:, None] >= -78) & (lat[:, None] <= -60)), ('ocean 45-60S', (lsm <= 0.5) & (lat[:, None] >= -60) & (lat[:, None] < -45)),
                ('all south of 60S', lat[:, None] < -60)):
    k = k & np.isfinite(G) & np.ones_like(P, bool); w = A[k]
    pm, gm = np.average(P[k], weights=w), np.average(G[k], weights=w)
    print(f'  {name:<18} model {pm:6.0f} mm/yr ({pm * w.sum() / 1e12 / 1000 * 1000 / 1e3 * 1e3 / 1e3:7.0f} Gt/yr)   GPCP {gm:6.0f} mm/yr   model/GPCP {pm / gm:.2f}')
