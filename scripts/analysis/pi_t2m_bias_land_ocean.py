"""Land/ocean split of the epoch-corrected PI T2m bias (companion of pi_t2m_bias_epoch_corrected.py).
ERA5 1990-2014 (annual 0.25 deg; seasonal 2 deg yseasmean) interpolated bilinearly to the model's remapped grid.
dT_obs = HadCRUT5 1990-2014 minus 1850-1900 per 5-deg box (>= 80 % of 1850-1900 months present), boxes classed
land/ocean by the model lsm averaged onto them (> 0.5 land). Usage: <outdata/oifs> <y0> <y1> <label>
"""
import sys, glob, numpy as np, xarray as xr
d, y0, y1, label = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
E = '/albedo/work/user/jstreffi/obs/era5/netcdf'
HC = '/albedo/work/user/jstreffi/obs/HadCRUT5/HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc'
def yrs(pat): return sorted(f for f in glob.glob(f'{d}/{pat}_*.nc') if y0 <= int(f[-12:-8]) <= y1)
m = xr.open_mfdataset(yrs('atm_remapped_1m_2t'), combine='by_coords', use_cftime=True)['2t']
td = 'time_counter' if 'time_counter' in m.dims else 'time'
lsm = xr.open_dataset(yrs('atm_remapped_1m_lsm')[0], use_cftime=True)['lsm'].isel({td: 0}).load()
la, lo = m.lat, m.lon
def onmodel(da):
    da = da.rename({k: v for k, v in {'latitude': 'lat', 'longitude': 'lon'}.items() if k in da.dims}).sortby('lat')
    if float(da.lon.min()) < 0: da = da.assign_coords(lon=da.lon % 360).sortby('lon')
    return da.interp(lat=la, lon=lo % 360 if float(lo.min()) < 0 else lo)
era_ann = xr.open_dataset(f'{E}/T2M_yearmean.nc'); v = [k for k in era_ann.data_vars if era_ann[k].ndim >= 3][0]
ERA = {'ANN': onmodel(era_ann[v].mean(era_ann[v].dims[0]))}
for s in ('DJF', 'JJA'): ERA[s] = onmodel(xr.open_dataset(f'{E}/T2M_{s}.nc')['tas'].squeeze(drop=True))
hc = xr.open_dataset(HC)['tas_mean']; hm, hy = hc.time.dt.month, hc.time.dt.year
lsm5 = lsm.assign_coords(lon=lsm.lon % 360).sortby('lon').interp(lat=hc.latitude, lon=hc.longitude % 360, method='nearest')
# nearest is a crude class for a 5-deg box; refine with the coarsened mean where the grids allow
SEAS = {'ANN': range(1, 13), 'DJF': [12, 1, 2], 'JJA': [6, 7, 8]}
BANDS = [('30N-45N', 30, 45), ('45N-60N', 45, 60), ('60N-75N', 60, 75), ('75N-90N', 75, 90), ('30S-30N', -30, 30), ('60S-30S', -60, -30)]
w = np.cos(np.deg2rad(la)) * xr.ones_like(lsm)
wh = np.cos(np.deg2rad(hc.latitude)) * xr.ones_like(lsm5)
print(f'{label} {y0}-{y1}: (model - ERA5 1990-2014) + dT_obs(HadCRUT5 1850-1900 -> 1990-2014), land = lsm > 0.5')
print(f'{"seas":4s} {"band":8s} {"sfc":5s} {"mod-ERA5":>9s} {"dT_obs":>7s} {"nbox":>5s} {"PI bias":>8s}')
for s, mo in SEAS.items():
    mm = m.sel({td: m[td].dt.month.isin(list(mo))}).mean(td).load() - 273.15
    diff = mm - (ERA[s] - 273.15)
    hs = hc.where(hm.isin(list(mo)), drop=True); hsy = hy.where(hm.isin(list(mo)), drop=True)
    early = hs.where((hsy >= 1850) & (hsy <= 1900), drop=True); late = hs.where((hsy >= 1990) & (hsy <= 2014), drop=True)
    dT = (late.mean('time') - early.mean('time')).where(early.notnull().mean('time') >= 0.8)
    for b, a, z in BANDS:
        for sfc, msk, msk5 in (('land', lsm > 0.5, lsm5 > 0.5), ('ocean', lsm <= 0.5, lsm5 <= 0.5)):
            sel = msk & (la >= a) & (la < z)
            db = float((diff * w).where(sel).sum() / w.where(sel & diff.notnull()).sum())
            s5 = msk5 & (hc.latitude >= a) & (hc.latitude < z) & dT.notnull()
            n5 = int(s5.sum()); dt = float((dT * wh).where(s5).sum() / wh.where(s5).sum()) if n5 else np.nan
            print(f'{s:4s} {b:8s} {sfc:5s} {db:+9.2f} {dt:+7.2f} {n5:5d} {db + dt:+8.2f}')
