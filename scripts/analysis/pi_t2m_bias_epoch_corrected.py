"""T2m bias of a PI-control run against a pre-industrial reference built from ERA5 + HadCRUT5.

  PI reference = ERA5(1990-2014) - dT_obs,  dT_obs = HadCRUT5(1990-2014) - HadCRUT5(1850-1900)
  PI bias      = model - PI reference = (model - ERA5) + dT_obs

dT_obs is taken per 5-degree HadCRUT5 box (analysis ensemble mean), then averaged over each band with
cos(lat) weights, over the boxes that have data in at least 80 % of the 1850-1900 months in that season
(coverage is reported; the polar bands rest on few boxes before 1900). Model and ERA5 band means use all
grid points, so dT_obs is a band mean applied to a full-band bias.
Usage: python pi_t2m_bias_epoch_corrected.py <outdata/oifs> <y0> <y1> <label>
"""
import sys, glob, numpy as np, xarray as xr
d, y0, y1, label = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
ERA5Y = '/albedo/work/user/jstreffi/obs/era5/netcdf/T2M_yearmean.nc'     # 1990-2014 annual
HC = '/albedo/work/user/jstreffi/obs/HadCRUT5/HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc'
BANDS = [('global', -90, 90), ('90S-60S', -90, -60), ('60S-30S', -60, -30), ('30S-30N', -30, 30),
         ('30N-60N', 30, 60), ('60N-90N', 60, 90), ('70N-90N', 70, 90)]
SEAS = {'ANN': list(range(1, 13)), 'DJF': [12, 1, 2], 'JJA': [6, 7, 8]}

def latname(da):
    return 'lat' if 'lat' in da.dims else 'latitude'

def band_mean(da, lo, hi):
    la = latname(da)
    sel = da.where((da[la] >= lo) & (da[la] <= hi))
    w = np.cos(np.deg2rad(da[la])) * xr.ones_like(sel)
    w = w.where(sel.notnull())
    return float((sel * w).sum() / w.sum())

files = sorted(f for f in glob.glob(f'{d}/atm_remapped_1m_2t_*.nc') if y0 <= int(f[-12:-8]) <= y1)
m = xr.open_mfdataset(files, combine='by_coords', use_cftime=True)['2t']
tdim = 'time_counter' if 'time_counter' in m.dims else 'time'
mon = m[tdim].dt.month

era = xr.open_dataset(ERA5Y)
ev = [v for v in era.data_vars if era[v].ndim >= 3][0]
era_ann = era[ev].mean(era[ev].dims[0])

hc = xr.open_dataset(HC)['tas_mean']
hmon = hc['time'].dt.month
hyr = hc['time'].dt.year

print(f'{label}: model {y0}-{y1} vs ERA5 1990-2014, PI reference via HadCRUT5 1850-1900 -> 1990-2014')
print(f'{"season":5s} {"band":9s} {"model":>8s} {"ERA5":>8s} {"mod-ERA5":>9s} {"dT_obs":>7s} {"cover":>6s} {"PI bias":>8s}')
for s, months in SEAS.items():
    mm = m.sel({tdim: mon.isin(months)}).mean(tdim).compute() - 273.15
    if s == 'ANN':
        ee = era_ann - 273.15
    else:   # yseasmean of ERA5 monthly 1990-2014, r180x91 (T2M_r.nc split by season)
        ee = xr.open_dataset(f'/albedo/work/user/jstreffi/obs/era5/netcdf/T2M_{s}.nc')['tas'].squeeze() - 273.15
    hs = hc.where(hmon.isin(months), drop=True)
    hsy = hyr.where(hmon.isin(months), drop=True)
    early = hs.where((hsy >= 1850) & (hsy <= 1900), drop=True)
    late = hs.where((hsy >= 1990) & (hsy <= 2014), drop=True)
    cover = early.notnull().mean('time')
    dT = (late.mean('time') - early.mean('time')).where(cover >= 0.8)
    for name, lo, hi in BANDS:
        mv, ev_ = band_mean(mm, lo, hi), band_mean(ee, lo, hi)
        dt = band_mean(dT, lo, hi)
        la = latname(cover)
        cb = cover.where((cover[la] >= lo) & (cover[la] <= hi))
        cov = float((cb >= 0.8).sum() / cb.notnull().sum())
        print(f'{s:5s} {name:9s} {mv:8.2f} {ev_:8.2f} {mv-ev_:+9.2f} {dt:+7.2f} {cov:6.0%} {mv-ev_+dt:+8.2f}')
