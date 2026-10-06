"""Two checks on the Southern Ocean surface cap of PICAL_crunveg_tke_albsn082.

1. Westerlies: zonal-mean 10 m zonal wind over the ocean, 30-75S, model 2195-2209 against ERA5
   (period of obs/era5/netcdf/u10.nc), annual and by season: jet maximum, its latitude, and the mean
   over 50-65S, where the wind-driven upwelling of the seasonal ice zone is set. ERA5 is present day;
   the jet has strengthened and moved poleward since pre-industrial (ozone, greenhouse gases), so a
   pre-industrial reference would be weaker and further north than ERA5.

2. Freshwater delivered by land ice: FESOM 'runoff' (rivers plus the Antarctic snow accumulation that
   the runoff mapper returns as calving) and the melt in the ice-shelf cavities ('fw' at cavity nodes),
   integrated south of 60S in Gt/yr, with the split of the runoff by latitude band and the share that
   falls into the regions NEW and OLD of so_destabilisation.py.

Usage: so_westerlies_and_freshwater.py   (reval environment; writes data/clim/so_westerlies_and_freshwater.txt)
"""
import os, glob
import numpy as np
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
RUN = f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082/outdata'
ARM = f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_kpp_gm1000/outdata/fesom'
Y0, Y1 = 2195, 2209
SEAS = {'ANN': list(range(1, 13)), 'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8], 'SON': [9, 10, 11]}
out = []

# ---- 1. westerlies -------------------------------------------------------------------------
fs = sorted(f for f in glob.glob(f'{RUN}/oifs/atm_remapped_1m_10u_*.nc') if Y0 <= int(f[-12:-8]) <= Y1)
m = xr.open_mfdataset(fs, combine='by_coords', use_cftime=True)['10u']
t = 'time_counter' if 'time_counter' in m.dims else 'time'
mc = m.groupby(f'{t}.month').mean(t).compute().sortby('lat')
lsm = xr.open_dataset(f'{RUN}/oifs/atm_remapped_1m_lsm_{Y1}-{Y1}.nc')['lsm'].isel({d: 0 for d in ('time', 'time_counter') if d in xr.open_dataset(f'{RUN}/oifs/atm_remapped_1m_lsm_{Y1}-{Y1}.nc')['lsm'].dims}).sortby('lat')
era = xr.open_dataset('/albedo/work/user/jstreffi/obs/era5/netcdf/u10.nc')
ev = [v for v in era.data_vars if era[v].ndim >= 3][0]
e = era[ev].rename({k: v for k, v in (('latitude', 'lat'), ('longitude', 'lon')) if k in era[ev].dims})
et = [d for d in e.dims if d not in ('lat', 'lon')][0]
eper = f"{str(e[et].values[0])[:7]} to {str(e[et].values[-1])[:7]}"
ec = e.groupby(f'{et}.month').mean(et).sortby('lat').interp(lat=mc['lat'], lon=mc['lon']).compute()
ocean = lsm < 0.5
out += [f'1. Zonal-mean 10 m zonal wind over the ocean, model {Y0}-{Y1} against ERA5 ({eper}).',
        f"{'season':7s}{'jet max model':>14s}{'ERA5':>8s}{'jet lat model':>15s}{'ERA5':>8s}{'mean 50-65S model':>19s}{'ERA5':>8s}{'ratio u^2':>11s}"]
for s, mm in SEAS.items():
    r = []
    for x in (mc, ec):
        z = x.sel(month=mm).mean('month').where(ocean).mean('lon').sel(lat=slice(-75, -30))
        zi = z.interp(lat=np.arange(-75, -30, 0.1))
        r.append((float(zi.max()), float(zi['lat'][int(zi.argmax())]), float(z.sel(lat=slice(-65, -50)).mean())))
    out.append(f"{s:7s}{r[0][0]:14.2f}{r[1][0]:8.2f}{r[0][1]:15.1f}{r[1][1]:8.1f}{r[0][2]:19.2f}{r[1][2]:8.2f}{(r[0][2] / r[1][2]) ** 2:11.2f}")

# ---- 2. freshwater from land ice ------------------------------------------------------------
md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lon, lat = md['lon'].values, md['lat'].values
area = md['nod_area'].values[0]
cav = md['ulevels_nod2D'].values > 1
od = (~cav) & (np.abs(md['zbar_n_bottom'].values) > 2000) & (lat < -55)
mld = lambda d: np.abs(xr.open_dataset(f'{d}/MLD2.fesom.2177.nc')['MLD2'].isel(time=8).values)
ma, mt = mld(ARM), mld(f'{RUN}/fesom')
REG = {'NEW': od & (ma > 1000) & (mt < 500), 'OLD': od & (ma > 1000) & (mt > 1000)}
GT = 1000.0 * 86400 * 365 / 1e12                      # m3/s of water -> Gt/yr
ro = np.zeros(len(lat)); fwc = np.zeros(len(lat)); n = 0
for y in range(2200, 2210):
    ro += np.nan_to_num(xr.open_dataset(f'{RUN}/fesom/runoff.fesom.{y}.nc')['runoff'].mean('time').values)
    fwc += np.nan_to_num(xr.open_dataset(f'{RUN}/fesom/fw.fesom.{y}.nc')['fw'].mean('time').values); n += 1
ro /= n; fwc /= n
S = lat < -60
out += ['', '2. Freshwater from land ice south of 60S, 2200-2209 mean [Gt/yr].',
        f"   runoff field (calving of the Antarctic snow accumulation + rivers): {np.sum((ro * area)[S]) * GT:8.0f}",
        f"   melt in the ice-shelf cavities (fw at cavity nodes):               {np.sum((fwc * area)[S & cav]) * GT:8.0f}",
        '   runoff by latitude band [Gt/yr] and share of the total south of 50S:']
tot = np.sum((ro * area)[lat < -50]) * GT
for a, b in ((-90, -75), (-75, -70), (-70, -65), (-65, -60), (-60, -55), (-55, -50)):
    v = np.sum((ro * area)[(lat >= a) & (lat < b)]) * GT
    out.append(f"     {a:4d} to {b:4d}: {v:8.0f}  {100 * v / tot:5.1f} %")
for k, r in REG.items():
    v = np.sum((ro * area)[r]) * GT
    out.append(f"   runoff into region {k} ({area[r].sum() / 1e12:.2f}e6 km2): {v:6.0f} Gt/yr = {1000 * v * 1e12 / 1000 / area[r].sum():.0f} mm/yr of freshwater")
nz = (ro > 0) & (lat < -50)
out.append(f"   nodes receiving runoff south of 50S: {nz.sum()} of {(lat < -50).sum()}; area {area[nz].sum() / 1e12:.2f}e6 km2; "
           f"90 % of it enters over {np.sort((ro * area)[nz])[::-1].cumsum().searchsorted(0.9 * (ro * area)[nz].sum()) + 1} nodes")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/so_westerlies_and_freshwater.txt', 'w').write(txt + '\n')
