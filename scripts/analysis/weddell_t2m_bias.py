"""T2m bias against ERA5 (1990-2014) in the Weddell Sea, for the two candidate lines.

The reval T2m bias maps show a warm patch over the Weddell Sea in both PICAL_crunveg_tke_albsn082
(TKE only) and PICAL_crunveg_idemix_albsn082 (TKE + IDEMIX).  This gives its size: mean over the
Weddell box (60W-0E, 78S-60S, ocean and ice shelf alike) and over its warmest part, the peak
grid-point value, and the area warmer than 3 K and 5 K, for the annual mean and the seasons.
Model 2195-2209, remapped monthly 2t; ERA5 interpolated to the model's regular grid.
ERA5 is present-day; a pre-industrial reference would be colder, so the PI bias is larger than
these numbers by the observed warming there (not applied: HadCRUT5 has almost no pre-1900 data here).

Usage: weddell_t2m_bias.py    (reval environment; writes data/clim/weddell_t2m_bias_2195-2209.txt)
"""
import os, glob
import numpy as np
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4'
E = '/albedo/work/user/jstreffi/obs/era5/netcdf/T2M.nc'
RUNS = [('TKE only', 'PICAL_crunveg_tke_albsn082'), ('TKE + IDEMIX', 'PICAL_crunveg_idemix_albsn082')]
Y0, Y1 = 2195, 2209
SEAS = {'ANN': list(range(1, 13)), 'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8], 'SON': [9, 10, 11]}
BOX = dict(lon0=-60, lon1=0, lat0=-78, lat1=-60)

era = xr.open_dataset(E)
ev = [v for v in era.data_vars if era[v].ndim >= 3][0]
et = era[ev].rename({'latitude': 'lat', 'longitude': 'lon'}) if 'latitude' in era[ev].dims else era[ev]
tname = [d for d in et.dims if d not in ('lat', 'lon')][0]
eclim = et.groupby(f'{tname}.month').mean(tname).sortby('lat')


def model_clim(run):
    fs = sorted(f for f in glob.glob(f'{P}/{run}/outdata/oifs/atm_remapped_1m_2t_*.nc') if Y0 <= int(f[-12:-8]) <= Y1)
    m = xr.open_mfdataset(fs, combine='by_coords', use_cftime=True)['2t']
    t = 'time_counter' if 'time_counter' in m.dims else 'time'
    return m.groupby(f'{t}.month').mean(t).compute().sortby('lat')


out = [f'T2m minus ERA5 (1990-2014) in the Weddell Sea, {BOX["lon0"]}..{BOX["lon1"]}E, {BOX["lat0"]}..{BOX["lat1"]}N; model {Y0}-{Y1}.',
       'box = area mean; core = mean over the warmest 25 % of the box area; peak = warmest grid point (with its position);',
       '>3K, >5K = share of the box area warmer than that.', '']
for label, run in RUNS:
    mc = model_clim(run)
    ec = eclim.interp(lat=mc['lat'], lon=mc['lon'] % 360 if float(eclim['lon'].max()) > 180 else mc['lon'])
    ec = ec.assign_coords(lon=mc['lon'])
    b = mc - ec
    lon180 = ((b['lon'] + 180) % 360) - 180
    b = b.assign_coords(lon=lon180).sortby('lon').sel(lon=slice(BOX['lon0'], BOX['lon1']), lat=slice(BOX['lat0'], BOX['lat1']))
    w = (np.cos(np.deg2rad(b['lat'])) * xr.ones_like(b.isel(month=0))).values
    out.append(f'{label}  ({run})')
    out.append(f"{'season':6s} {'box':>6s} {'core':>6s} {'peak':>6s} {'at':>14s} {'>3K':>6s} {'>5K':>6s}")
    for s, mm in SEAS.items():
        x = b.sel(month=mm).mean('month')
        v = x.values
        ok = np.isfinite(v)
        box = np.sum(v[ok] * w[ok]) / np.sum(w[ok])
        order = np.argsort(v[ok])[::-1]
        cw = np.cumsum(w[ok][order]); top = order[cw <= 0.25 * cw[-1]]
        core = np.sum(v[ok][top] * w[ok][top]) / np.sum(w[ok][top])
        i, j = np.unravel_index(np.nanargmax(v), v.shape)
        f3 = np.sum(w[ok & (v > 3)]) / np.sum(w[ok]); f5 = np.sum(w[ok & (v > 5)]) / np.sum(w[ok])
        out.append(f"{s:6s} {box:6.2f} {core:6.2f} {v[i, j]:6.2f} {float(x['lat'][i]):6.1f}N {float(x['lon'][j]):6.1f}E {100 * f3:5.0f}% {100 * f5:5.0f}%")
    out.append('')
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/weddell_t2m_bias_{Y0}-{Y1}.txt', 'w').write(txt + '\n')
