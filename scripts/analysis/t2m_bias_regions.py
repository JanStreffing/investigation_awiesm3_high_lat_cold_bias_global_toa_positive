"""Where is a coupled run warmest against ERA5?  T2m bias by band, surface and season, plus a map.

T2m (OpenIFS remapped monthly 2t) over Y0-Y1 minus the ERA5 monthly climatology (T2M.nc),
on the model's remapped grid.  Bias per latitude band x {land, ocean} x {ANN, DJF, JJA},
cos-lat weighted, and a two-panel annual/DJF-JJA map.  ERA5 is an IFS sibling but is fine
for T2m (campaign protocol).  For a pre-industrial run the target is a COLD bias against
present-day ERA5: roughly the observed 1850 -> 1990-2020 warming, about -1 K globally with
larger values at high northern latitudes.

Usage:  ARM=PI200 Y0=1480 Y1=1499 python3 scripts/analysis/t2m_bias_regions.py
        ROOT=/path/to/runtime ... (run outside the default root)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); Y0, Y1 = int(os.environ.get('Y0', 1480)), int(os.environ.get('Y1', 1499))
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'
BANDS = [('90-60N', 60, 90), ('60-45N', 45, 60), ('45-30N', 30, 45), ('30N-30S', -30, 30),
         ('30-45S', -45, -30), ('45-60S', -60, -45), ('60-90S', -90, -60)]


def load(var, y):
    p = f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(p):
        p = f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
    with xr.open_dataset(p, decode_times=False) as d:
        k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values).astype('f8'), np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)


t2 = []
for y in range(Y0, Y1 + 1):
    a, lat, lon = load('2t', y); t2.append(a)
T = np.stack(t2).mean(0)                                   # (12, lat, lon)
T = T - 273.15 if np.nanmean(T) > 100 else T
with xr.open_dataset(LSMF, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values); m = m[0] if m.ndim == 3 else m
land = m > 0.5
with xr.open_dataset('/work/ab0246/a270092/obs/era5/netcdf/T2M.nc') as d:
    c = d['t2m'].groupby('time.month').mean('time')
    la = [x for x in c.dims if 'lat' in x][0]; lo = [x for x in c.dims if 'lon' in x][0]
    O = np.asarray(c.interp({la: ('lat', lat), lo: ('lon', np.sort(lon % 360))}).values, float)
O = O - 273.15 if np.nanmean(O) > 100 else O
order = np.argsort(lon % 360); T = T[:, :, order]; land = land[:, order]
L = np.broadcast_to(lat[:, None], land.shape); W = np.cos(np.deg2rad(L))
B = {'ANN': (T - O).mean(0), 'DJF': (T - O)[[11, 0, 1]].mean(0), 'JJA': (T - O)[[5, 6, 7]].mean(0)}
am = lambda f, k: float(np.average(f[k & np.isfinite(f)], weights=W[k & np.isfinite(f)]))
print(f'{ARM} {Y0}-{Y1}: T2m minus ERA5 [K]')
print(f'  global ANN {am(B["ANN"], np.ones_like(land)):+.2f}   land {am(B["ANN"], land):+.2f}   ocean {am(B["ANN"], ~land):+.2f}')
print(f'  {"band":<9}' + ''.join(f'{s + " " + x:>12}' for s in ('ANN', 'DJF', 'JJA') for x in ('land', 'ocean')))
for name, lo_, hi_ in BANDS:
    k = (L >= lo_) & (L < hi_)
    print(f'  {name:<9}' + ''.join(f'{am(B[s], k & (land if x == "land" else ~land)):12.2f}' for s in ('ANN', 'DJF', 'JJA') for x in ('land', 'ocean')))

LON = np.sort(lon % 360)
fig, axes = plt.subplots(3, 1, figsize=(9, 11), dpi=150)
for ax, s in zip(axes, ('ANN', 'DJF', 'JJA')):
    im = ax.pcolormesh(LON, lat, B[s], cmap='RdBu_r', vmin=-6, vmax=6, shading='auto')
    ax.contour(LON, lat, land.astype(float), levels=[0.5], colors='k', linewidths=0.4)
    ax.set_title(f'({"abc"[("ANN","DJF","JJA").index(s)]}) {s} T2m minus ERA5, global {am(B[s], np.ones_like(land)):+.2f} K', loc='left', fontsize=10)
    ax.set_ylabel('lat')
fig.colorbar(im, ax=axes, fraction=0.025, label='K')
fig.suptitle(f'{ARM} {Y0}-{Y1} (pre-industrial target: about -1 K against present-day ERA5)', x=0.01, ha='left')
out = os.path.join(REPO, 'plots', f't2m_bias_regions_{ARM}_{Y0}-{Y1}.png'); fig.savefig(out); print('saved', out)
