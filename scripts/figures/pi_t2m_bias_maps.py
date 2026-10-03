"""Maps and zonal means of the epoch-corrected PI T2m bias (see scripts/analysis/pi_t2m_bias_epoch_corrected.py).
  PI bias = (model - ERA5 1990-2014) + dT_obs,  dT_obs = HadCRUT5 1990-2014 minus 1850-1900 (analysis ensemble mean)
dT_obs is used per 5-deg box where >= 80 % of the 1850-1900 months exist; elsewhere it is filled with the mean of
the covered boxes in the same latitude row (stippled), and with the global covered mean where a row has none.
Usage: python pi_t2m_bias_maps.py <outdata/oifs> <y0> <y1> <label> <out.png>
"""
import sys, glob, numpy as np, xarray as xr, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm
d, y0, y1, label, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
E = '/albedo/work/user/jstreffi/obs/era5/netcdf'
HC = '/albedo/work/user/jstreffi/obs/HadCRUT5/HadCRUT.5.0.2.0.analysis.anomalies.ensemble_mean.nc'
def yrs(pat): return sorted(f for f in glob.glob(f'{d}/{pat}_*.nc') if y0 <= int(f[-12:-8]) <= y1)
m = xr.open_mfdataset(yrs('atm_remapped_1m_2t'), combine='by_coords', use_cftime=True)['2t']
td = 'time_counter' if 'time_counter' in m.dims else 'time'
lsm = xr.open_dataset(yrs('atm_remapped_1m_lsm')[0], use_cftime=True)['lsm'].isel({td: 0}).load()
la, lo = m.lat, m.lon
lo360 = lo % 360

def onmodel(da, method='linear'):
    da = da.rename({k: v for k, v in {'latitude': 'lat', 'longitude': 'lon'}.items() if k in da.dims}).sortby('lat')
    da = da.assign_coords(lon=da.lon % 360).sortby('lon')
    da = xr.concat([da.isel(lon=-1).assign_coords(lon=da.lon[-1] - 360), da, da.isel(lon=0).assign_coords(lon=da.lon[0] + 360)], 'lon')
    return da.interp(lat=la, lon=lo360, method=method, kwargs={'fill_value': None})

era_ann = xr.open_dataset(f'{E}/T2M_yearmean.nc'); v = [k for k in era_ann.data_vars if era_ann[k].ndim >= 3][0]
ERA = {'ANN': onmodel(era_ann[v].mean(era_ann[v].dims[0]))}
for s in ('DJF', 'JJA'): ERA[s] = onmodel(xr.open_dataset(f'{E}/T2M_{s}.nc')['tas'].squeeze(drop=True))
hc = xr.open_dataset(HC)['tas_mean']; hm, hy = hc.time.dt.month, hc.time.dt.year
SEAS = {'ANN': range(1, 13), 'DJF': [12, 1, 2], 'JJA': [6, 7, 8]}
w = np.cos(np.deg2rad(la))

lev = np.array([-6, -4, -3, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 3, 4, 6])
cmap = plt.get_cmap('RdBu_r', len(lev) + 1); norm = BoundaryNorm(lev, cmap.N, extend='both')
fig = plt.figure(figsize=(15, 13))
gs = fig.add_gridspec(3, 2, width_ratios=[3.2, 1], wspace=0.08, hspace=0.18)
for i, (s, mo) in enumerate(SEAS.items()):
    mm = m.sel({td: m[td].dt.month.isin(list(mo))}).mean(td).load() - 273.15
    raw = mm - (ERA[s] - 273.15)
    hs = hc.where(hm.isin(list(mo)), drop=True); hsy = hy.where(hm.isin(list(mo)), drop=True)
    early = hs.where((hsy >= 1850) & (hsy <= 1900), drop=True); late = hs.where((hsy >= 1990) & (hsy <= 2014), drop=True)
    covered = early.notnull().mean('time') >= 0.8
    dT = (late.mean('time') - early.mean('time')).where(covered)
    wl = np.cos(np.deg2rad(dT.latitude)) * xr.ones_like(dT)
    gmean = float((dT * wl).sum() / wl.where(dT.notnull()).sum())
    dTf = dT.fillna(dT.mean('longitude')).fillna(gmean)
    dTm = onmodel(dTf, 'nearest')
    filled = onmodel((~covered).astype(float), 'nearest') > 0.5
    pib = raw + dTm
    glob_pib = float((pib * w).sum() / (w * xr.ones_like(pib)).sum())
    glob_raw = float((raw * w).sum() / (w * xr.ones_like(raw)).sum())

    ax = fig.add_subplot(gs[i, 0], projection=ccrs.Robinson())
    def cyc(a):   # wrap column by hand: the remapped lon axis is not exactly equally spaced
        return np.concatenate([a, a[:, :1]], 1)
    lonc = np.concatenate([lo360.values, lo360.values[:1] + 360]); z = cyc(pib.values)
    cf = ax.contourf(lonc, la, z, levels=lev, cmap=cmap, norm=norm, extend='both', transform=ccrs.PlateCarree())
    fz = cyc(filled.astype(float).values)
    ax.contourf(lonc, la, fz, levels=[0.5, 1.5], colors='none', hatches=['...'], transform=ccrs.PlateCarree())
    ax.coastlines(lw=0.4, color='#333333'); ax.set_global()
    ax.set_title(f'{s}: PI bias = (model - ERA5 1990-2014) + HadCRUT5 warming     global {glob_pib:+.2f} K '
                 f'(vs ERA5 {glob_raw:+.2f})', fontsize=10, loc='left')

    az = fig.add_subplot(gs[i, 1])
    def zm(x, sel=None):
        x = x if sel is None else x.where(sel)
        return x.mean('lon')
    az.axvline(0, color='0.6', lw=0.8)
    az.plot(zm(raw), la, color='0.45', lw=1.2, ls='--', label='model - ERA5')
    az.plot(zm(dTm), la, color='tab:orange', lw=1.2, label='HadCRUT5 warming')
    az.plot(zm(pib, lsm > 0.5), la, color='tab:green', lw=1.6, label='PI bias, land')
    az.plot(zm(pib, lsm <= 0.5), la, color='tab:blue', lw=1.6, label='PI bias, ocean')
    az.plot(zm(pib), la, color='k', lw=2.0, label='PI bias, all')
    az.set_ylim(-90, 90); az.set_yticks(range(-90, 91, 30)); az.set_xlim(-5, 5); az.grid(alpha=0.3)
    az.yaxis.tick_right(); az.set_xlabel('K', fontsize=9); az.set_title(f'{s} zonal mean', fontsize=10, loc='left')
    if i == 0: az.legend(fontsize=7.5, loc='lower left')
cb = fig.colorbar(cf, ax=fig.axes[0::2], orientation='horizontal', fraction=0.025, pad=0.03, ticks=lev)
cb.set_label(f'T2m bias against a pre-industrial reference, K  ({label}, {y0}-{y1})\n'
             'stippled: no HadCRUT5 1850-1900 coverage, latitude-row mean warming used', fontsize=9)
fig.savefig(out, dpi=120, bbox_inches='tight')
print('written', out)
