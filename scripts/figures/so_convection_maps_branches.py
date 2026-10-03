"""Where the Southern Ocean convects, where the winter warm bias sits, and where a branch acts.

Southern Ocean (50-80S) maps over Y0-Y1: JJA mixed-layer depth (MLD2) of CONTROL and each
branch; CONTROL's JJA T2m bias against ERA5; each branch's JJA T2m change; JJA effective
ice thickness of CONTROL and the branch change; and the 200-500 m temperature of CONTROL
against PHC3 (PHC3 binned from the mesh nodes).  If a lever removes the convection but the
bias sits elsewhere, this is where it shows.

Usage:  CONTROL=PI200 RUNS=PI200_spp,PI200_h0 Y0=1565 Y1=1579 \
        python3 scripts/figures/so_convection_maps_branches.py
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
MESH = '/work/ab0246/a270092/input/fesom2/core3/fesom.mesh.diag.nc'
PHC = '/work/ab0246/a270092/postprocessing/climatologies/CORE3_220509'
JJA = [5, 6, 7]


def gr(arm, v, months):
    acc = 0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{arm}/outdata/fesom/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
            a = d[v].values.astype('f8')
            if a.shape[0] > 12:
                n = [31, 29 if a.shape[0] == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; e = np.cumsum([0] + n)
                a = np.stack([np.nanmean(a[e[i]:e[i + 1]], 0) for i in range(12)])
            acc = acc + a[months].mean(0) / (Y1 - Y0 + 1)
            la = d['lat'].values; lo = d['lon'].values
    acc[np.abs(acc) > 1e30] = np.nan
    return acc, la, lo


def t2m(arm):
    acc = 0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{arm}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc', decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            acc = acc + np.squeeze(d[k].values)[JJA].mean(0).astype('f8') / (Y1 - Y0 + 1)
            la = np.squeeze(d['lat'].values); lo = np.squeeze(d['lon'].values)
    o = np.argsort(lo % 360)
    return acc[:, o], la, (lo % 360)[o]


T, alat, alon = t2m(CONTROL)
with xr.open_dataset('/work/ab0246/a270092/obs/era5/netcdf/T2M.nc') as d:
    c = d['t2m'].groupby('time.month').mean('time')
    la = [x for x in c.dims if 'lat' in x][0]; lo = [x for x in c.dims if 'lon' in x][0]
    c = c.assign_coords({lo: c[lo] % 360}).sortby(lo)
    E = np.asarray(c.interp({la: ('lat', alat), lo: ('lon', alon)}).values, float)[JJA].mean(0)
MLD = {a: gr(a, 'MLD2', JJA)[0] for a in [CONTROL] + RUNS}
HI = {a: gr(a, 'm_ice', JJA)[0] for a in [CONTROL] + RUNS}
_, glat, glon = gr(CONTROL, 'MLD2', JJA)
T300, _, _ = gr(CONTROL, 'temp', list(range(12)))
with xr.open_dataset(f'{R}/{CONTROL}/outdata/fesom/temp.fesom.gr.{Y0}.nc', decode_times=False) as d:
    z = d['nz'].values
kz = (z >= 200) & (z < 500)
T300 = np.nanmean(T300[..., kz], -1)
with xr.open_dataset(MESH) as m:
    nlon = m['lon'].values % 360; nlat = m['lat'].values; area = m['nod_area'].values[0]
with xr.open_dataset(f'{PHC}/temp.fesom.1958.nc', decode_times=False) as d:
    p = np.squeeze(d['temp'].values)[:, :47].astype('f8'); p[np.abs(p) > 1e30] = np.nan
p300 = np.nanmean(p[:, kz], 1)
iy = np.clip(((nlat + 90) / 0.5).astype(int), 0, len(glat) - 1); ix = np.clip((nlon / 0.5).astype(int), 0, len(glon) - 1)
ok = np.isfinite(p300)
num = np.zeros((len(glat), len(glon))); den = np.zeros_like(num)
np.add.at(num, (iy[ok], ix[ok]), p300[ok] * area[ok]); np.add.at(den, (iy[ok], ix[ok]), area[ok])
P300 = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)

n = len(RUNS)
fig, axs = plt.subplots(3 + n, 2, figsize=(14, 3.0 * (3 + n)), constrained_layout=True)
def panel(ax, lon, lat, f, vmin, vmax, cmap, title):
    im = ax.pcolormesh(lon, lat, f, vmin=vmin, vmax=vmax, cmap=cmap, shading='auto'); ax.set_ylim(-80, -50); ax.set_title(title, fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.8); return im
panel(axs[0, 0], glon, glat, -MLD[CONTROL], 0, 800, 'viridis', f'(a) {CONTROL} JJA mixed layer depth [m]')
panel(axs[0, 1], alon, alat, T - E, -6, 6, 'RdBu_r', f'(b) {CONTROL} JJA T2m minus ERA5 [K]')
panel(axs[1, 0], glon, glat, HI[CONTROL], 0, 2, 'Blues', f'(c) {CONTROL} JJA effective ice thickness [m]')
panel(axs[1, 1], glon, glat, T300 - P300, -2, 2, 'RdBu_r', f'(d) {CONTROL} 200-500 m temperature minus PHC3 [K]')
panel(axs[2, 0], glon, glat, P300, -2, 3, 'RdYlBu_r', '(e) PHC3 200-500 m temperature [C]')
panel(axs[2, 1], glon, glat, T300, -2, 3, 'RdYlBu_r', f'(f) {CONTROL} 200-500 m temperature [C]')
for i, a in enumerate(RUNS):
    Ta, _, _ = t2m(a)
    panel(axs[3 + i, 0], glon, glat, -MLD[a], 0, 800, 'viridis', f'({chr(103 + 2 * i)}) {a} JJA mixed layer depth [m]')
    panel(axs[3 + i, 1], alon, alat, Ta - T, -3, 3, 'RdBu_r', f'({chr(104 + 2 * i)}) {a} minus {CONTROL}, JJA T2m [K]')
fig.suptitle(f'Southern Ocean, {Y0}-{Y1}: convection, winter warm bias, ice thickness, deep reservoir, and branch responses')
out = os.path.join(REPO, 'plots', f'so_convection_maps_{Y0}-{Y1}.png'); fig.savefig(out, dpi=100); print('saved', out)
