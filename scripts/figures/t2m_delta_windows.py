"""T2m change between two windows of one run, with a no-change reference, as maps and band means.

Panel (a): mean(A0-A1) minus mean(B0-B1), the change being tested.
Panel (b): mean(B0-B1) minus mean(C0-C1), the preceding window-to-window difference with no
namelist change, so drift and decadal variability can be read off against (a).
Annual means of OpenIFS remapped monthly 2t.  Band means are cos-lat weighted, land and ocean
separately.  Default: PI200, RCL_INPPMIN 70000 -> 50000 (S4) at 1500.

Usage:  python3 scripts/figures/t2m_delta_windows.py
        ARM=PI200 A=1500:1509 B=1490:1499 C=1480:1489 LABEL="S4 on" python3 ...
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); LABEL = os.environ.get('LABEL', 'RCL_INPPMIN 70000 -> 50000 (S4)')
win = lambda k, d: tuple(int(x) for x in os.environ.get(k, d).split(':'))
A, B, C = win('A', '1500:1509'), win('B', '1490:1499'), win('C', '1480:1489')
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
BANDS = [('90-60N', 60, 90), ('60-45N', 45, 60), ('45-30N', 30, 45), ('30N-30S', -30, 30),
         ('30-45S', -45, -30), ('45-60S', -60, -45), ('60-90S', -90, -60)]


def mean_t2(y0, y1):
    acc = []
    for y in range(y0, y1 + 1):
        with xr.open_dataset(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc', decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            acc.append(np.squeeze(d[k].values).astype('f8').mean(0))
            lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
    return np.mean(acc, 0), lat, lon


TA, lat, lon = mean_t2(*A); TB, _, _ = mean_t2(*B); TC, _, _ = mean_t2(*C)
with xr.open_dataset(LSMF, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values); m = m[0] if m.ndim == 3 else m
o = np.argsort(lon % 360); LON = (lon % 360)[o]
D1, D0, land = (TA - TB)[:, o], (TB - TC)[:, o], (m > 0.5)[:, o]
L = np.broadcast_to(lat[:, None], land.shape); W = np.cos(np.deg2rad(L))
am = lambda f, k: float(np.average(f[k], weights=W[k]))
g = lambda f: am(f, np.ones_like(land))
print(f'{ARM}: (a) {A[0]}-{A[1]} minus {B[0]}-{B[1]} [{LABEL}]   (b) {B[0]}-{B[1]} minus {C[0]}-{C[1]} [no change]')
print(f'  global: (a) {g(D1):+.3f} K   (b) {g(D0):+.3f} K')
print(f'  {"band":<9}{"(a) land":>10}{"(a) ocean":>11}{"(b) land":>10}{"(b) ocean":>11}')
for nm, lo_, hi_ in BANDS:
    k = (L >= lo_) & (L < hi_)
    print(f'  {nm:<9}{am(D1, k & land):10.2f}{am(D1, k & ~land):11.2f}{am(D0, k & land):10.2f}{am(D0, k & ~land):11.2f}')

fig, axes = plt.subplots(2, 1, figsize=(9, 8.5), dpi=150)
for ax, D, t in ((axes[0], D1, f'(a) {A[0]}-{A[1]} minus {B[0]}-{B[1]}: {LABEL}, global {g(D1):+.2f} K'),
                 (axes[1], D0, f'(b) {B[0]}-{B[1]} minus {C[0]}-{C[1]}: no change, global {g(D0):+.2f} K')):
    im = ax.pcolormesh(LON, lat, D, cmap='RdBu_r', vmin=-2, vmax=2, shading='auto')
    ax.contour(LON, lat, land.astype(float), levels=[0.5], colors='k', linewidths=0.4)
    ax.set_title(t, loc='left', fontsize=10); ax.set_ylabel('lat')
fig.colorbar(im, ax=axes, fraction=0.03, label='annual-mean T2m change [K]')
fig.suptitle(f'{ARM}: decade-mean T2m change', x=0.01, ha='left')
out = os.path.join(REPO, 'plots', f't2m_delta_{ARM}_{A[0]}-{A[1]}_vs_{B[0]}-{B[1]}.png'); fig.savefig(out); print('saved', out)
