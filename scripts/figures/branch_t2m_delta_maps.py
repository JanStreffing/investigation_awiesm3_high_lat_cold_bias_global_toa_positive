"""T2m change of each branch against its control over the same years: annual and JJA maps.

Branch minus CONTROL, OpenIFS remapped monthly 2t, mean over Y0-Y1, one row per branch,
columns ANN and JJA.  Cells are stippled where the difference is smaller than twice the
standard error of the difference of two Y0-Y1 means estimated from the control's
year-to-year standard deviation at that cell (no autocorrelation correction, so the
unstippled area is a generous reading).  Titles carry the global and 60-90S means.

Usage:  CONTROL=PI200 RUNS=PI200_spp,PI200_mle,PI200_h0,PI200_all3 Y0=1570 Y1=1579 \
        python3 scripts/figures/branch_t2m_delta_maps.py
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
Y0, Y1 = int(os.environ.get('Y0', 1570)), int(os.environ.get('Y1', 1579))
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'


def load(arm):
    ann, jja = [], []
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{arm}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc', decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8'); lat = np.squeeze(d['lat'].values); lon = np.squeeze(d['lon'].values)
        ann.append(a.mean(0)); jja.append(a[[5, 6, 7]].mean(0))
    return np.array(ann), np.array(jja), lat, lon


C_ann, C_jja, lat, lon = load(CONTROL)
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = lsm[0] if lsm.ndim == 3 else lsm
W = np.cos(np.deg2rad(lat))[:, None] * np.ones(len(lon))[None, :]
N = Y1 - Y0 + 1
fig, axs = plt.subplots(len(RUNS), 2, figsize=(13, 2.9 * len(RUNS)), constrained_layout=True)
for i, arm in enumerate(RUNS):
    B_ann, B_jja, _, _ = load(arm)
    for j, (lab, b, c) in enumerate((('ANN', B_ann, C_ann), ('JJA', B_jja, C_jja))):
        dlt = b.mean(0) - c.mean(0)
        se = 2 * c.std(0, ddof=1) * np.sqrt(2.0 / N)
        ax = axs[i, j]
        im = ax.pcolormesh(lon, lat, dlt, cmap='RdBu_r', vmin=-2 if lab == 'ANN' else -3, vmax=2 if lab == 'ANN' else 3, shading='auto')
        ax.contour(lon, lat, lsm, levels=[0.5], colors='k', linewidths=0.4)
        yy, xx = np.meshgrid(np.arange(0, len(lat), 4), np.arange(0, len(lon), 4), indexing='ij')
        m = np.abs(dlt[yy, xx]) < se[yy, xx]
        ax.scatter(lon[xx][m], lat[yy][m], s=0.15, c='0.45', lw=0)
        g = np.average(dlt, weights=W); sp = np.average(dlt[lat < -60], weights=W[lat < -60])
        ax.set_title(f'({chr(97 + 2 * i + j)}) {arm} - {CONTROL}, {lab} {Y0}-{Y1}:  global {g:+.2f} K, 60-90S {sp:+.2f} K', fontsize=9)
        ax.set_ylim(-90, 90)
        if j == 1 or i == len(RUNS) - 1:
            pass
    fig.colorbar(im, ax=axs[i, :], shrink=0.9, label='K')
fig.suptitle(f'T2m change against {CONTROL}, {Y0}-{Y1} (stippled: within 2 standard errors of the control spread)')
out = os.path.join(REPO, 'plots', f'branch_t2m_delta_{Y0}-{Y1}.png'); fig.savefig(out, dpi=110); print('saved', out)
