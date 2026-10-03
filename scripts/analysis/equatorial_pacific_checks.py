"""Cheap checks on the equatorial Pacific cold tongue: winds, convection position, SST gradient.

Along 5S-5N, per longitude band, model last-decade means against observations:
  u10    10 m zonal wind [m/s]          model 10u      vs ERA5 u10.nc 'uas', 1989-2014 (dynamics: ERA5 is fine)
  precip [mm/day]                       model tp       vs GPCP precip.mon.mean.nc (1980-2009)
  OLR    [W/m2]                         model -ttr     vs CERES EBAF toa_lw_all (2000s)
  SST    [degC]                         model sst      vs HadISST 1980-2009
Plus the eastern edge of the convecting warm pool: the easternmost longitude (130E-280E)
where precip > 6 mm/day, and where OLR < 240 W/m2.
IFS accumulations (tp in m, ttr in J/m2) are per hourly output step.
Only 1990-forced arms are compared with present-day observations.

Usage:  python3 scripts/analysis/equatorial_pacific_checks.py
        RUNS=16E_1990:1390:1399,11V:1380:1389 python3 ...
        RUNS=/abs/path/to/run:1650:1679 ...   (a run outside ROOT)
        FIG=plots/equatorial_pacific_profiles.png ...   (also draw the profiles)
"""
import os
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
OBS = '/work/ab0246/a270092/obs'
RUNS = [r.split(':') for r in os.environ.get('RUNS', '16E_1990:1390:1399,11V:1380:1389').split(',')]
BANDS = [('west 130-160E', 130, 160), ('Nino4 160E-150W', 160, 210), ('Nino3.4 170-120W', 190, 240),
         ('Nino3 150-90W', 210, 270)]


def model_mean(run, var, y0, y1):
    acc = None
    base = run if run.startswith('/') else f'{R}/{run}'
    for y in range(int(y0), int(y1) + 1):
        p = f'{base}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(p):   # some runs write a doubled frequency infix
            p = f'{base}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        with xr.open_dataset(p, decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8').mean(0)
            lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values) % 360
        acc = a if acc is None else acc + a
    return acc / (int(y1) - int(y0) + 1), lat, lon


def eq_profile(f, lat, lon):
    """5S-5N cos-weighted mean as a function of longitude (sorted 0-360)."""
    k = (lat >= -5) & (lat <= 5)
    w = np.cos(np.deg2rad(lat[k]))
    prof = np.nansum(f[k] * w[:, None], 0) / np.sum(w)
    o = np.argsort(lon)
    return lon[o], prof[o]


def band(lonp, prof, a, b):
    k = (lonp >= a) & (lonp <= b)
    return float(np.nanmean(prof[k]))


def east_edge(lonp, prof, test):
    k = (lonp >= 130) & (lonp <= 280) & test(prof)
    return float(lonp[k].max()) if k.any() else np.nan


def obs_field(path, var, t0=None, t1=None, scale=1.0, flip_lon=False):
    with xr.open_dataset(path) as d:
        v = d[var]
        if t0 and 'time' in v.dims:
            v = v.sel(time=slice(t0, t1))
        a = v.mean([c for c in v.dims if 'time' in c]).values.astype('f8') * scale
        latn = [c for c in d.coords if c.lower().startswith('lat')][0]
        lonn = [c for c in d.coords if c.lower().startswith('lon')][0]
        return a, d[latn].values, d[lonn].values % 360


rows = {}
for run, y0, y1 in RUNS:
    u, la, lo = model_mean(run, '10u', y0, y1)
    tp, _, _ = model_mean(run, 'tp', y0, y1)
    ttr, _, _ = model_mean(run, 'ttr', y0, y1)
    sst, _, _ = model_mean(run, 'sst', y0, y1)
    sst = np.where(sst > 200, sst - 273.15, np.nan)
    rows[f'{os.path.basename(run)} {y0}-{y1}'] = {'u10': eq_profile(u, la, lo), 'precip': eq_profile(tp * 1000 * 24, la, lo),
                                'OLR': eq_profile(-ttr / 3600, la, lo), 'SST': eq_profile(sst, la, lo)}

u, la, lo = obs_field(f'{OBS}/era5/netcdf/u10.nc', 'uas')   # ERA5 10 m zonal wind, 1989-2014
pr, pla, plo = obs_field(f'{OBS}/gpcp/precip.mon.mean.nc', 'precip', '1980-01-01', '2009-12-31')
with xr.open_dataset(f'{OBS}/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc', decode_times=False) as c:
    olr, cla, clo = c['toa_lw_all_clim'].values.mean(0), c['lat'].values, c['lon'].values % 360
with xr.open_dataset(f'{OBS}/hadisst2/HadISST_sst.nc') as h:
    s = h['sst'].sel(time=slice('1980-01-01', '2009-12-31'))
    hs = s.where(s > -100).mean('time').values; hla, hlo = h['latitude'].values, h['longitude'].values % 360
rows['OBS (ERA5/GPCP/CERES/HadISST)'] = {'u10': eq_profile(u, la, lo), 'precip': eq_profile(pr, pla, plo),
                                         'OLR': eq_profile(olr, cla, clo), 'SST': eq_profile(hs, hla, hlo)}

for var, unit in (('u10', 'm/s'), ('precip', 'mm/day'), ('OLR', 'W/m2'), ('SST', 'degC')):
    print(f'\n=== {var} [{unit}], 5S-5N')
    print(f'{"":30s}' + ''.join(f'{b[0]:>18s}' for b in BANDS) + ('   east edge' if var in ('precip', 'OLR') else ''))
    for name, d in rows.items():
        lonp, prof = d[var]
        line = f'{name:30s}' + ''.join(f'{band(lonp, prof, a, b):18.2f}' for _, a, b in BANDS)
        if var == 'precip':
            line += f'   {east_edge(lonp, prof, lambda p: p > 6):6.0f}E (>6 mm/d)'
        if var == 'OLR':
            line += f'   {east_edge(lonp, prof, lambda p: p < 240):6.0f}E (<240)'
        print(line)
for name, d in rows.items():
    lonp, prof = d['SST']
    print(f'SST gradient west(130-160E) minus Nino3 (150-90W), {name}: {band(lonp, prof, 130, 160) - band(lonp, prof, 210, 270):.2f} K')

FIG = os.environ.get('FIG', '')
if FIG:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    SURF, TXT1, TXT2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
    colours = ['#2a78d6', '#eb6834', '#1baf7a', '#e87ba4']
    fig, axes = plt.subplots(4, 1, figsize=(8.5, 10), dpi=150, sharex=True); fig.patch.set_facecolor(SURF)
    for ax, (var, unit) in zip(axes, (('u10', 'm s$^{-1}$'), ('precip', 'mm day$^{-1}$'), ('OLR', 'W m$^{-2}$'), ('SST', '$^\\circ$C'))):
        ax.set_facecolor(SURF)
        for i, (name, d) in enumerate(rows.items()):
            lonp, prof = d[var]; k = (lonp >= 120) & (lonp <= 280)
            obs = name.startswith('OBS')
            ax.plot(lonp[k], prof[k], color='k' if obs else colours[i % len(colours)], lw=2.2 if obs else 1.6,
                    ls='--' if obs else '-', label='observations' if obs else name)
        ax.set_ylabel(f'{var} [{unit}]', color=TXT2, fontsize=9); ax.grid(color=GRID, lw=0.8)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        ax.tick_params(colors=TXT2, labelsize=8.5)
    axes[0].legend(fontsize=8, frameon=False, loc='lower left')
    axes[-1].set_xlabel('longitude [deg E], mean over 5S-5N', color=TXT2, fontsize=9)
    fig.suptitle('Equatorial Pacific along the equator: 10 m zonal wind, rain, OLR, SST', x=0.01, ha='left', fontsize=11, color=TXT1)
    fig.tight_layout(); fig.savefig(FIG, facecolor=SURF); print('saved', FIG)
