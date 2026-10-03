"""Month-by-month surface state over the Southern Ocean pack: when does the warm bias form?

Over cells where both the model (ci) and HadISST (sic) have >= 0.8 ice in that month,
55-80S ocean, on the 1x1 grid: T2m and skin bias against ERA5 (1990-2014), the model's
inversion (T2m minus skin) against ERA5's, clear-sky and all-sky LW down against CERES,
total cloud cover, sensible heat flux, and the pack area.  IFS fluxes accumulated J/m2
per hour (/3600), positive down.

Usage:  ARM=PI200 Y0=1565 Y1=1579 python3 scripts/analysis/so_pack_seasonal_cycle.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); Y0, Y1 = int(os.environ.get('Y0', 1565)), int(os.environ.get('Y1', 1579))
CER = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
SIC = '/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc'
E5 = '/work/ab0246/a270092/obs/era5/netcdf'
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'

with xr.open_dataset(CER, decode_times=False) as dc:
    clat = np.squeeze(dc['lat'].values); clon = np.squeeze(dc['lon'].values)
    C = {k: np.squeeze(dc[v].values) for k, v in (('LWdn', 'sfc_lw_down_all_clim'), ('LWdn_clr', 'sfc_lw_down_clr_c_clim'))}
with xr.open_dataset(SIC, decode_times=False) as ds:
    sic = np.squeeze(ds['sic'].values); slat = np.squeeze(ds['latitude'].values)
if np.nanmax(sic) > 1.5: sic = sic / 100.0
if slat[0] > slat[-1]: sic = sic[:, ::-1, :]
sic = np.nan_to_num(sic, nan=0.0)


def era5(f, v):
    with xr.open_dataset(f'{E5}/{f}') as d:
        c = d[v].sel(time=d['time.year'].isin(range(1990, 2015))).groupby('time.month').mean('time')
        c = c.assign_coords(longitude=c['longitude'] % 360).sortby('longitude').sortby('latitude')
        return np.asarray(c.interp(latitude=('y', clat), longitude=('x', clon)).values, float)


O = {'2t': era5('T2M.nc', 't2m'), 'skt': era5('SKT_mon.nc', 'skt')}
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = lsm[0] if lsm.ndim == 3 else lsm
    mlat = np.squeeze(d['lat'].values); mlon = np.squeeze(d['lon'].values)
iy = np.abs(clat[:, None] - mlat[None, :]).argmin(1); ix = np.abs(((clon[:, None] - mlon[None, :] + 180) % 360 - 180)).argmin(1)
rg = lambda a: a[..., iy, :][..., ix]
LAT = np.broadcast_to(clat[:, None], sic.shape[1:]); W = np.cos(np.deg2rad(LAT))
OCN = (rg(lsm) <= 0.5) & (LAT >= -80) & (LAT <= -55)
VARS = ('2t', 'skt', 'ci', 'strd', 'strdc', 'sshf', 'tcc'); ACC = {'strd', 'strdc', 'sshf'}
M = {v: 0.0 for v in VARS}
for y in range(Y0, Y1 + 1):
    for v in VARS:
        with xr.open_dataset(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc', decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8')
        M[v] = M[v] + rg(a / 3600.0 if v in ACC else a) / (Y1 - Y0 + 1)
print(f'{ARM} {Y0}-{Y1}, 55-80S, cells with model ci >= 0.8 and HadISST >= 0.8 in that month')
print(f'  {"mon":>3}{"Mkm2":>6}{"T2m bias":>9}{"skin bias":>10}{"inv mod":>8}{"inv ERA5":>9}{"LWclr-CER":>10}{"LWall-CER":>10}{"tcc":>6}{"sshf":>6}{"ci mod":>7}')
for m in range(12):
    ci = np.nan_to_num(M['ci'][m]); k = OCN & (ci >= 0.8) & (sic[m] >= 0.8) & np.isfinite(M['2t'][m])
    if k.sum() < 10:
        print(f'  {m + 1:3d}  (fewer than 10 cells)'); continue
    av = lambda f: float(np.average(f[k], weights=W[k]))
    print(f'  {m + 1:3d}{float((W * 111.195e3 ** 2)[k].sum() / 1e12):6.2f}{av(M["2t"][m] - O["2t"][m]):+9.2f}{av(M["skt"][m] - O["skt"][m]):+10.2f}'
          f'{av(M["2t"][m] - M["skt"][m]):+8.2f}{av(O["2t"][m] - O["skt"][m]):+9.2f}{av(M["strdc"][m] - C["LWdn_clr"][m]):+10.1f}'
          f'{av(M["strd"][m] - C["LWdn"][m]):+10.1f}{av(M["tcc"][m]):6.2f}{av(M["sshf"][m]):6.1f}{av(ci):7.2f}')
