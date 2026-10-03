"""Are the Southern Hemisphere westerlies too strong or too far south?  Model against ERA5.

Zonal-mean 10 m zonal wind (model 10u, ERA5 uas_vas.nc 1990-2014) and zonal-mean zonal
wind on pressure levels (model pl_u, ERA5 U_mon.nc climatology), annual and JJA/DJF, over
Y0-Y1.  Reports the surface jet maximum and its latitude, the mean wind in 50-65S (the
Ekman-upwelling band for CDW) and the 850/500/200 hPa jet.  Stronger or more poleward
westerlies mean more CDW upwelled toward the continent.

Usage:  ARM=PI200 Y0=1565 Y1=1579 python3 scripts/analysis/so_westerlies_vs_era5.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); Y0, Y1 = int(os.environ.get('Y0', 1565)), int(os.environ.get('Y1', 1579))
E5 = '/work/ab0246/a270092/obs/era5/netcdf'
SEAS = {'ANN': list(range(12)), 'JJA': [5, 6, 7], 'DJF': [11, 0, 1]}


def zm(var):
    acc = 0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc', decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            acc = acc + np.squeeze(d[k].values).astype('f8').mean(-1) / (Y1 - Y0 + 1)
            lat = np.squeeze(d['lat'].values); plev = np.squeeze(d['pressure_levels'].values) if 'pressure_levels' in d else None
    return acc, lat, plev


U10, lat, _ = zm('10u'); UP, _, plev = zm('pl_u')
with xr.open_dataset(f'{E5}/uas_vas.nc') as d:
    e = d['u10'].sel(time=d['time.year'].isin(range(1990, 2015))).groupby('time.month').mean('time').mean('longitude').sortby('latitude')
    EU10 = np.asarray(e.interp(latitude=lat).values, float)
with xr.open_dataset(f'{E5}/U_mon.nc') as d:
    e = d['U'].mean('lon').sortby('lat')
    EUP = np.asarray(e.interp(lat=lat, plev=plev).values, float)
W = np.cos(np.deg2rad(lat)); band = lambda f, lo, hi: float(np.average(f[(lat >= lo) & (lat <= hi)], weights=W[(lat >= lo) & (lat <= hi)]))
sh = lat < -30
for s, ms in SEAS.items():
    m, o = U10[ms].mean(0), EU10[ms].mean(0)
    print(f'\n{s} {ARM} {Y0}-{Y1} against ERA5 1990-2014, zonal mean')
    print(f'  10 m zonal wind: jet max model {m[sh].max():.2f} m/s at {lat[sh][m[sh].argmax()]:.1f}, ERA5 {o[sh].max():.2f} at {lat[sh][o[sh].argmax()]:.1f}')
    for lo, hi in ((-65, -50), (-70, -60), (-50, -40)):
        print(f'    {lo}..{hi}: model {band(m, lo, hi):6.2f}  ERA5 {band(o, lo, hi):6.2f}  diff {band(m, lo, hi) - band(o, lo, hi):+.2f}')
    for p in (85000, 50000, 20000):
        i = int(np.argmin(np.abs(plev - p))); mm, oo = UP[ms].mean(0)[i], EUP[ms].mean(0)[i]
        print(f'  {p / 100:5.0f} hPa: jet max model {mm[sh].max():.1f} at {lat[sh][mm[sh].argmax()]:.1f}, ERA5 {oo[sh].max():.1f} at {lat[sh][oo[sh].argmax()]:.1f};'
              f'  50-65S diff {band(mm, -65, -50) - band(oo, -65, -50):+.2f} m/s')
