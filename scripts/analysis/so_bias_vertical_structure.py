"""Is the Southern Ocean winter warm bias a boundary-layer feature or a tropospheric one?

Zonal-mean temperature on pressure levels, CONTROL (Y0-Y1) minus the ERA5 monthly
climatology (T_mon.nc, 19 levels), for JJA and DJF, plus each branch minus CONTROL.  A
surface-driven bias (ocean heat through thin ice) is confined below about 850 hPa and
capped by the inversion; an atmosphere-driven one (heat transport, free-troposphere
warmth) reaches 500 hPa.  Levels below the surface are extrapolated by both products and
are read with that in mind (over the ice zone the surface is at about 990 hPa).

Prints band means at each level for 60-90S, 65-75S and 60-90N, writes a two-row figure:
bias JJA / DJF (top) and spp change JJA / DJF (bottom).

Usage:  CONTROL=PI200 RUNS=PI200_spp,PI200_h0 Y0=1565 Y1=1579 \
        python3 scripts/analysis/so_bias_vertical_structure.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
CONTROL = os.environ.get('CONTROL', 'PI200')
RUNS = os.environ.get('RUNS', 'PI200_spp,PI200_h0').split(',')
Y0, Y1 = int(os.environ.get('Y0', 1565)), int(os.environ.get('Y1', 1579))
SEAS = {'JJA': [5, 6, 7], 'DJF': [11, 0, 1]}
SHOW = [100000, 92500, 85000, 70000, 50000, 30000]


def load(arm):
    acc = 0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{arm}/outdata/oifs/atm_remapped_1m_pl_t_{y}-{y}.nc', decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8')                    # (12, plev, lat, lon)
            acc = acc + a.mean(-1) / (Y1 - Y0 + 1)
            lat = np.squeeze(d['lat'].values); plev = np.squeeze(d['pressure_levels'].values)
    return acc, lat, plev                                                 # (12, plev, lat)


C, lat, plev = load(CONTROL)
with xr.open_dataset('/work/ab0246/a270092/obs/era5/netcdf/T_mon.nc') as d:
    e = d['t'] if 't' in d else d[[v for v in d.data_vars if v not in ('time_bnds',)][0]]
    e = e.mean('lon').sortby('lat')
    E = np.asarray(e.interp(lat=lat, plev=plev).values, float)           # (12, plev, lat)
    if np.nanmean(E) < 100: E = E + 273.15
W = np.cos(np.deg2rad(lat))
B = {a: load(a)[0] for a in RUNS}
band = lambda f, lo, hi: float(np.average(f[(lat >= lo) & (lat <= hi)], weights=W[(lat >= lo) & (lat <= hi)]))

fig, axs = plt.subplots(1 + len(RUNS), 2, figsize=(12, 3.2 * (1 + len(RUNS))), constrained_layout=True)
for j, (s, ms) in enumerate(SEAS.items()):
    bias = C[ms].mean(0) - E[ms].mean(0)
    print(f'\n{s} {CONTROL} {Y0}-{Y1} minus ERA5, zonal mean [K]      then branch minus {CONTROL}')
    print(f'  {"hPa":>6}{"60-90S":>9}{"65-75S":>9}{"45-60S":>9}{"60-90N":>9}   ' + ''.join(f'{a.replace(CONTROL + "_", ""):>24}' for a in RUNS))
    print(f'  {"":>42}   ' + ''.join(f'{"60-90S":>8}{"65-75S":>8}{"60-90N":>8}' for a in RUNS))
    for p in SHOW:
        i = int(np.argmin(np.abs(plev - p)))
        row = f'  {p / 100:6.0f}{band(bias[i], -90, -60):9.2f}{band(bias[i], -75, -65):9.2f}{band(bias[i], -60, -45):9.2f}{band(bias[i], 60, 90):9.2f}   '
        for a in RUNS:
            dl = B[a][ms].mean(0) - C[ms].mean(0)
            row += f'{band(dl[i], -90, -60):8.2f}{band(dl[i], -75, -65):8.2f}{band(dl[i], 60, 90):8.2f}'
        print(row)
    ax = axs[0, j]
    im = ax.contourf(lat, plev / 100, bias, levels=np.arange(-5, 5.1, 0.5), cmap='RdBu_r', extend='both')
    ax.set_yscale('log'); ax.invert_yaxis(); ax.set_ylim(1000, 100); ax.set_title(f'({chr(97 + j)}) {CONTROL} minus ERA5, {s}, zonal mean T [K]', fontsize=9)
    ax.set_ylabel('hPa')
    for i, a in enumerate(RUNS):
        ax = axs[1 + i, j]
        dl = B[a][ms].mean(0) - C[ms].mean(0)
        im2 = ax.contourf(lat, plev / 100, dl, levels=np.arange(-1.5, 1.51, 0.15), cmap='RdBu_r', extend='both')
        ax.set_yscale('log'); ax.invert_yaxis(); ax.set_ylim(1000, 100); ax.set_ylabel('hPa')
        ax.set_title(f'({chr(97 + 2 * (i + 1) + j)}) {a} minus {CONTROL}, {s} [K]', fontsize=9)
fig.colorbar(im, ax=axs[0, :], shrink=0.8, label='bias [K]'); fig.colorbar(im2, ax=axs[1:, :], shrink=0.5, label='change [K]')
out = os.path.join(REPO, 'plots', f'so_bias_vertical_structure_{Y0}-{Y1}.png'); fig.savefig(out, dpi=110); print('saved', out)
