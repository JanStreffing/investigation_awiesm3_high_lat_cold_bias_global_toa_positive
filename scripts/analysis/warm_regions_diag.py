"""What drives the two warm regions of the PI bias map (scripts/figures/pi_t2m_bias_maps.py)?

A. NH summer land (JJA): central Asia/steppe against Europe and North America. Model vs ERA5 T2m
   (1990-2014), GPCP precipitation (1990-2014), CERES EBAF surface SW down and cloud area (2000-2021).
   Plus the model's own partition: SSR, SLHF, SSHF, evaporative fraction, top-layer soil water, LAI.
   CERES cldarea is a MODIS mask, not a radiative cover, so only its sign is read against tcc.
B. Marginal ice zones: every ocean point classed by model ice (ci) and HadISST 2.2 (1990-2014) ice,
   threshold 0.15, and the model-ERA5 T2m bias averaged per class. SH in JJA, NH in DJF.
All biases are model minus present-day obs; the PI epoch correction (HadCRUT5) is not applied here.
OpenIFS accumulations are per output hour: J/m2 -> /3600 W/m2, m -> *24000 mm/day.
Usage: python warm_regions_diag.py <outdata/oifs> <y0> <y1> <label>
"""
import sys, glob, numpy as np, xarray as xr
d, y0, y1, label = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
O = '/albedo/work/user/jstreffi/obs'
def yrs(v): return sorted(f for f in glob.glob(f'{d}/atm_remapped_1m_{v}_*.nc') if y0 <= int(f[-12:-8]) <= y1)
def mod(v, months):
    da = xr.open_mfdataset(yrs(v), combine='by_coords', use_cftime=True)[v]
    td = 'time_counter' if 'time_counter' in da.dims else 'time'
    return da.sel({td: da[td].dt.month.isin(months)}).mean(td).load()
t2 = mod('2t', [6, 7, 8]); la, lo = t2.lat, t2.lon % 360
def onmodel(da, method='linear'):
    da = da.rename({k: v for k, v in {'latitude': 'lat', 'longitude': 'lon'}.items() if k in da.dims}).sortby('lat')
    da = da.assign_coords(lon=da.lon % 360).sortby('lon')
    da = xr.concat([da.isel(lon=-1).assign_coords(lon=da.lon[-1] - 360), da, da.isel(lon=0).assign_coords(lon=da.lon[0] + 360)], 'lon')
    return da.interp(lat=la, lon=lo, method=method, kwargs={'fill_value': None})
lsm = xr.open_dataset(yrs('lsm')[0], use_cftime=True)['lsm'].squeeze(drop=True)
w = np.cos(np.deg2rad(la)) * xr.ones_like(t2)
def bmean(x, sel): return float((x * w).where(sel).sum() / w.where(sel & x.notnull()).sum())

# ---------------- A. summer land ----------------
JJA = [6, 7, 8]
M = {'T2m': t2 - 273.15, 'precip': mod('tp', JJA) * 24000, 'SW down': mod('ssrd', JJA) / 3600,
     'tcc': mod('tcc', JJA) * 100, 'SSR': mod('ssr', JJA) / 3600, 'LH': -mod('slhf', JJA) / 3600,
     'SH': -mod('sshf', JJA) / 3600, 'swvl1': mod('swvl1', JJA), 'LAI': mod('lai_lv', JJA) + mod('lai_hv', JJA)}
M['EF'] = M['LH'] / (M['LH'] + M['SH'])
era = xr.open_dataset(f'{O}/era5/netcdf/T2M_JJA.nc')['tas'].squeeze(drop=True) - 273.15
gp = xr.open_dataset(f'{O}/gpcp/precip.mon.mean.nc')['precip']
gp = gp.sel(time=slice('1990', '2014')); gp = gp.where(gp.time.dt.month.isin(JJA), drop=True).mean('time')
ce = xr.open_dataset(f'{O}/CERES/CERES_EBAF_Ed4.1_Subset_200003-202106_JJA.nc').squeeze(drop=True)
OBS = {'T2m': onmodel(era), 'precip': onmodel(gp), 'SW down': onmodel(ce['sfc_sw_down_all_mon']),
       'tcc': onmodel(ce['cldarea_total_daynight_mon'])}
lon180 = ((lo + 180) % 360) - 180
BOX = {'C Asia 40-55N 50-110E': (40, 55, 50, 110), 'Europe 45-60N 0-40E': (45, 60, 0, 40),
       'N America 40-55N 115-85W': (40, 55, -115, -85), 'Siberia 55-75N 60-180E': (55, 75, 60, 180)}
print(f'{label} {y0}-{y1}  A. JJA land, model / obs / model-obs (obs: ERA5 T2m 90-14, GPCP 90-14, CERES 00-21)')
print(f'{"box":26s}' + ''.join(f'{k:>20s}' for k in ['T2m', 'precip', 'SW down', 'tcc']) +
      ''.join(f'{k:>8s}' for k in ['SSR', 'LH', 'SH', 'EF', 'swvl1', 'LAI']))
for b, (a, z, l0, l1) in BOX.items():
    sel = (lsm > 0.5) & (la >= a) & (la < z) & (lon180 >= l0) & (lon180 < l1)
    row = f'{b:26s}'
    for k in ['T2m', 'precip', 'SW down', 'tcc']:
        mv, ov = bmean(M[k], sel), bmean(OBS[k], sel)
        row += f'{mv:7.1f}/{ov:5.1f}/{mv-ov:+6.1f}'
    row += ''.join(f'{bmean(M[k], sel):8.2f}' for k in ['SSR', 'LH', 'SH', 'EF', 'swvl1', 'LAI'])
    print(row)

# ---------------- B. marginal ice zones ----------------
hs = xr.open_dataset(f'{O}/hadisst2/HadISST.2.2.0.0_sea_ice_concentration.nc', use_cftime=True)['sic']
hs = hs.sel(time=slice('1990', '2014'))
ERA = {'JJA': era + 0, 'DJF': xr.open_dataset(f'{O}/era5/netcdf/T2M_DJF.nc')['tas'].squeeze(drop=True) - 273.15}
print(f'\nB. marginal ice zones: model-ERA5 T2m by ice class (threshold 0.15; obs ice = HadISST 2.2 1990-2014)')
print(f'{"case":18s} {"class":26s} {"area %":>7s} {"mod-ERA5":>9s} {"ci mod":>7s} {"ci obs":>7s}')
for hem, s, months, a, z in (('SH', 'JJA', [6, 7, 8], -80, -50), ('NH', 'DJF', [12, 1, 2], 45, 90)):
    ci = mod('ci', months)
    tm = (t2 if s == 'JJA' else mod('2t', months)) - 273.15
    bias = tm - onmodel(ERA[s])
    ho = hs.where(hs.time.dt.month.isin(months), drop=True).mean('time')
    ho = onmodel(ho.where(ho >= 0), 'nearest')
    if float(ho.max()) > 1.5: ho = ho / 100
    ocean = (lsm <= 0.5) & (la >= a) & (la <= z) & ci.notnull() & ho.notnull()
    cls = {'both ice (pack)': (ci >= 0.15) & (ho >= 0.15), 'obs ice, model open': (ci < 0.15) & (ho >= 0.15),
           'model ice, obs open': (ci >= 0.15) & (ho < 0.15), 'both open': (ci < 0.15) & (ho < 0.15)}
    tot = float(w.where(ocean).sum())
    for c, msk in cls.items():
        sel = ocean & msk
        print(f'{hem+" "+s:18s} {c:26s} {100*float(w.where(sel).sum())/tot:7.1f} {bmean(bias, sel):+9.2f} '
              f'{bmean(ci, sel):7.2f} {bmean(ho, sel):7.2f}')
