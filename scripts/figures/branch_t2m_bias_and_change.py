"""T2m bias against present-day ERA5 for a control and its branches, and each branch's change.

One figure per season (ANN, DJF, JJA).  Rows: CONTROL, then each branch.  Left column: T2m
minus the ERA5 1990-2014 monthly climatology (on the model's remapped grid), mean over
Y0-Y1.  Right column: branch minus CONTROL over the same years, sparsely stippled where the
change is within 2*sd*sqrt(2/N) of the control's year-to-year spread at that cell (no
autocorrelation correction).  Titles carry the global, 60-90S and 60-90N means.  For a
pre-industrial run the expected bias against present-day ERA5 is negative, about -1 K
globally.

Usage:  CONTROL=PI200 RUNS=PI200_spp,PI200_mle,PI200_h0,PI200_all3 Y0=1565 Y1=1579 \
        python3 scripts/figures/branch_t2m_bias_and_change.py
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
RUNS = os.environ.get('RUNS', 'PI200_spp,PI200_mle,PI200_h0,PI200_all3').split(',')
Y0, Y1 = int(os.environ.get('Y0', 1565)), int(os.environ.get('Y1', 1579))
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'
SEAS = {'ANN': list(range(12)), 'DJF': [11, 0, 1], 'JJA': [5, 6, 7]}
LABEL = {'PI200_spp': 'spp on', 'PI200_mle': 'Fox-Kemper MLE', 'PI200_h0': 'h0min 0.25', 'PI200_all3': 'spp + MLE + h0min'}


def load(arm):
    out = []
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{arm}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc', decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            out.append(np.squeeze(d[k].values).astype('f8'))
            lat = np.squeeze(d['lat'].values); lon = np.squeeze(d['lon'].values)
    return np.array(out), lat, lon                                   # (years, 12, lat, lon)


C, lat, lon = load(CONTROL)
order = np.argsort(lon % 360); lon = (lon % 360)[order]; C = C[..., order]
D = {a: load(a)[0][..., order] for a in RUNS}
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = (lsm[0] if lsm.ndim == 3 else lsm)[:, order]
with xr.open_dataset('/work/ab0246/a270092/obs/era5/netcdf/T2M.nc') as d:
    c = d['t2m'].groupby('time.month').mean('time')
    la = [x for x in c.dims if 'lat' in x][0]; lo = [x for x in c.dims if 'lon' in x][0]
    c = c.assign_coords({lo: c[lo] % 360}).sortby(lo)
    E = np.asarray(c.interp({la: ('lat', lat), lo: ('lon', lon)}).values, float)      # (12, lat, lon)
W = np.cos(np.deg2rad(lat))[:, None] * np.ones(len(lon))[None, :]
N = Y1 - Y0 + 1
wm = lambda f, k: float(np.average(f[k & np.isfinite(f)], weights=W[k & np.isfinite(f)]))
ALL = np.ones(W.shape, bool); SP = (lat < -60)[:, None] & ALL; NP = (lat > 60)[:, None] & ALL
yy, xx = np.meshgrid(np.arange(0, len(lat), 6), np.arange(0, len(lon), 6), indexing='ij')

for s, ms in SEAS.items():
    cs = C[:, ms].mean(1)                                             # (years, lat, lon)
    obs = E[ms].mean(0)
    rows = [CONTROL] + RUNS
    fig, axs = plt.subplots(len(rows), 2, figsize=(14, 2.75 * len(rows)), constrained_layout=True)
    for i, arm in enumerate(rows):
        bs = cs if arm == CONTROL else D[arm][:, ms].mean(1)
        bias = bs.mean(0) - obs
        ax = axs[i, 0]
        im0 = ax.pcolormesh(lon, lat, bias, cmap='RdBu_r', vmin=-6, vmax=6, shading='auto')
        ax.contour(lon, lat, lsm, levels=[0.5], colors='k', linewidths=0.4)
        name = arm if arm == CONTROL else f'{arm} ({LABEL.get(arm, "")})'
        ax.set_title(f'({chr(97 + 2 * i)}) {name} minus ERA5, {s}:  global {wm(bias, ALL):+.2f}, '
                     f'60-90S {wm(bias, SP):+.2f}, 60-90N {wm(bias, NP):+.2f} K', fontsize=8.5)
        ax = axs[i, 1]
        if arm == CONTROL:
            ax.axis('off')
            ax.text(0.5, 0.5, f'{Y0}-{Y1} means ({N} years)\nleft: bias against ERA5 1990-2014\n'
                    f'right: branch minus {CONTROL}\nstippled: change within 2 sd sqrt(2/N)\nof the control spread',
                    ha='center', va='center', fontsize=10, transform=ax.transAxes)
            continue
        dlt = bs.mean(0) - cs.mean(0)
        thr = 2 * cs.std(0, ddof=1) * np.sqrt(2.0 / N)
        im1 = ax.pcolormesh(lon, lat, dlt, cmap='RdBu_r', vmin=-2, vmax=2, shading='auto')
        ax.contour(lon, lat, lsm, levels=[0.5], colors='k', linewidths=0.4)
        m = np.abs(dlt[yy, xx]) < thr[yy, xx]
        ax.scatter(lon[xx][m], lat[yy][m], s=0.3, c='0.35', lw=0)
        ax.set_title(f'({chr(97 + 2 * i + 1)}) {arm} minus {CONTROL}, {s}:  global {wm(dlt, ALL):+.2f}, '
                     f'60-90S {wm(dlt, SP):+.2f}, 60-90N {wm(dlt, NP):+.2f} K', fontsize=8.5)
    for ax in axs.flat:
        ax.set_ylim(-90, 90); ax.set_xticks(range(0, 361, 60)); ax.tick_params(labelsize=7)
    fig.colorbar(im0, ax=axs[:, 0], shrink=0.5, label='bias against ERA5 [K]', location='bottom')
    fig.colorbar(im1, ax=axs[:, 1], shrink=0.5, label=f'change against {CONTROL} [K]', location='bottom')
    fig.suptitle(f'T2m {s}, {Y0}-{Y1}: bias against present-day ERA5 (left) and change by each branch (right)')
    out = os.path.join(REPO, 'plots', f'branch_t2m_bias_and_change_{s}_{Y0}-{Y1}.png')
    fig.savefig(out, dpi=100); plt.close(fig); print('saved', out)
