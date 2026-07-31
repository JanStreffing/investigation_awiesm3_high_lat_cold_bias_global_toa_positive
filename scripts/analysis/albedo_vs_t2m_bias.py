"""Does the surface-albedo bias explain the land cold bias?

albedo_by_region.py showed the model is too bright over high-latitude land (Siberian tundra
+0.094, Canada/Siberia boreal ~+0.03, Sahara +0.020) while temperate and tropical vegetation
is within +-0.004. The question here is whether that error is large enough, and spatially
organised enough, to account for the land cold bias.

WEIGHTING. A raw annual-mean albedo difference is meaningless at high latitude, because
polar-winter albedo multiplies almost no sunlight. The physically meaningful quantity is the
ABSORBED SHORTWAVE the albedo error costs:

    dSW_alb(month) = -(alpha_model - alpha_CERES) * SWdown_model(month)

averaged over months. That is insolation-weighted by construction: polar winter contributes
~0 because SWdown ~ 0, and every tropical month counts fully. No seasonal masking needed.

PERIOD. amip_presentday (1990-2014) against ERA5 1990-2014 and the CERES 2005-2015
climatology, so the T2m reference carries no epoch offset (the control would need the
+0.42 K correction measured earlier).

Correlation is over land gridcells, area-weighted, 60S-75N (excluding permanent ice, whose
albedo and temperature are both dominated by the ice sheet rather than by vegetation).
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
OBS = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
ERA5 = '/work/ab0246/a270092/obs/era5/netcdf/T2M.nc'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
RUN, Y0, Y1 = 'amip_presentday', 1990, 2014
ACC = 3600.0

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values
_d = xr.open_dataset(f'{RT}/{RUN}/outdata/oifs/atm_remapped_1m_ssr_1m_{Y0}-{Y0}.nc')
lat, lon = _d['ssr'].lat.values, _d['ssr'].lon.values
_d.close()


def regrid(da):
    da = xr.concat([da.isel(lon=[-1]).assign_coords(lon=da.lon.values[-1:] - 360.), da,
                    da.isel(lon=[0]).assign_coords(lon=da.lon.values[:1] + 360.)], dim='lon')
    return da.interp(lat=xr.DataArray(np.clip(lat, -89.5, 89.5), dims='y'),
                     lon=xr.DataArray(np.where(lon < 0, lon + 360, lon), dims='x')).values


def ceres(v):
    return regrid(xr.open_dataset(OBS)[v])


def model(v):
    acc = []
    for y in range(Y0, Y1 + 1):
        d = xr.open_dataset(f'{RT}/{RUN}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc')
        a = d[v].values
        acc.append(a / ACC if v != '2t' else a)
        d.close()
    return np.mean(acc, axis=0)


# ERA5 monthly climatology on the model grid
_e = xr.open_dataset(ERA5)
_v = 't2m' if 't2m' in _e else list(_e.data_vars)[0]
_da = _e[_v].groupby('time.month').mean('time')
_da = _da.rename({('latitude' if 'latitude' in _da.dims else 'lat'): 'lat',
                  ('longitude' if 'longitude' in _da.dims else 'lon'): 'lon'}).sortby('lat')
e_t2m = np.stack([regrid(_da.isel(month=m)) for m in range(12)])
_e.close()

m_net, m_dn, m_t2m = model('ssr'), model('ssrd'), model('2t')
c_up, c_dn = ceres('sfc_sw_up_all_clim'), ceres('sfc_sw_down_all_clim')

a_m = (m_dn - m_net) / np.maximum(m_dn, 1e-6)
a_o = c_up / np.maximum(c_dn, 1e-6)
dSW = -np.mean((a_m - a_o) * m_dn, axis=0)          # insolation-weighted, W/m2
dT = np.mean(m_t2m - e_t2m, axis=0)                  # K, annual mean

land = (lsm > 0.5) & (lat[:, None] >= -60) & (lat[:, None] <= 75)
land &= np.isfinite(dSW) & np.isfinite(dT)
# drop Greenland: ice sheet, not vegetation
gl = (lat[:, None] >= 59) & (lat[:, None] <= 84) & \
     ((((lon + 180) % 360) - 180)[None, :] >= -73) & ((((lon + 180) % 360) - 180)[None, :] <= -12)
land &= ~gl

w = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], dSW.shape)[land]
x, y = dSW[land], dT[land]
r = np.corrcoef(x, y)[0, 1]
slope = np.polyfit(x, y, 1)[0]

print(f'{RUN} {Y0}-{Y1} vs ERA5/CERES, land 60S-75N excl. Greenland, {land.sum()} cells\n')
print(f'  area-weighted mean albedo-induced SW loss : {np.average(x, weights=w):+7.2f} W/m2')
print(f'  area-weighted mean T2m bias               : {np.average(y, weights=w):+7.2f} K')
print(f'  spatial correlation r                     : {r:+7.3f}')
print(f'  regression slope                          : {slope:+7.3f} K per W/m2\n')

BINS = [(-99, -8), (-8, -4), (-4, -2), (-2, -1), (-1, 0), (0, 99)]
print(f"  {'albedo-induced SW loss':26s}{'cells':>8s}{'mean dSW':>10s}{'mean dT2m':>11s}")
for lo, hi in BINS:
    m = (x >= lo) & (x < hi)
    if m.sum() < 20:
        continue
    print(f'  {f"{lo:+.0f} to {hi:+.0f} W/m2":26s}{m.sum():8d}'
          f'{np.average(x[m], weights=w[m]):10.2f}{np.average(y[m], weights=w[m]):11.2f}')
print('\n  dSW < 0 = model too bright, losing absorbed SW.  If the albedo error drives the')
print('  cold bias, cells losing more SW should be colder (positive slope).')
