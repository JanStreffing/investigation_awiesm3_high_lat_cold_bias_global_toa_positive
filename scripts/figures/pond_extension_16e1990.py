"""16E_1990 with the stronger melt ponds from 1400: sea ice and the pond-albedo chain, year by year.

Panels: (a) NH September and SH February sea-ice extent [M km2] from OpenIFS ci (cells with
ci > 0.15, cos-lat area on the remapped grid), (b) July FESOM melt-pond fraction and July
OpenIFS surface albedo over Arctic pack ice (a_ice / ci > 0.8, north of 70N).
The pond settings changed at 1400 (albpnd 0.28 -> 0.35, rfracmax 0.75 -> 0.6, pndaspect 1.3 -> 1.6).

Usage:  python3 scripts/figures/pond_extension_16e1990.py      (Y0=1370 Y1=1419)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = '/work/bb1469/a270092/runtime/awiesm3-v3.4/16E_1990/outdata'
Y0, Y1 = int(os.environ.get('Y0', 1370)), int(os.environ.get('Y1', 1419))
with xr.open_dataset('/work/ab0246/a270092/input/fesom2/core3/mesh.nc') as m:
    AREA = m['cell_area'].values.astype('f8')


def oifs(var, y):
    with xr.open_dataset(f'{R}/oifs/atm_remapped_1m_{var}_{y}-{y}.nc', decode_times=False) as d:
        k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values), np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)


def fesom_july(var, y):
    with xr.open_dataset(f'{R}/fesom/{var}.fesom.{y}.nc', decode_times=False) as d:
        v = np.squeeze(d[var].values); lat = d['lat'].values
    if v.shape[0] == 12:
        return v[6], lat
    j0 = 181 + (1 if v.shape[0] == 366 else 0)
    return v[j0:j0 + 31].mean(0), lat


years, nhsep, shfeb, apnd, fal = [], [], [], [], []
for y in range(Y0, Y1 + 1):
    ci, lat, lon = oifs('ci', y)
    L = np.broadcast_to(lat[:, None], ci.shape[1:])
    cell = 6.371e6 ** 2 * np.cos(np.deg2rad(L)) * np.deg2rad(abs(lat[1] - lat[0])) * 2 * np.pi / ci.shape[2]
    ext = lambda c, hemi: float(np.nansum(np.where((c > 0.15) & hemi, cell, 0.0))) / 1e12
    nhsep.append(ext(ci[8], L > 0)); shfeb.append(ext(ci[1], L < 0))
    fa, _, _ = oifs('fal', y)
    kk = (ci[6] > 0.8) & (L > 70) & np.isfinite(fa[6])
    fal.append(np.average(fa[6][kk], weights=np.cos(np.deg2rad(L[kk]))))
    a, flat = fesom_july('a_ice', y); p, _ = fesom_july('apnd', y)
    k = (a > 0.8) & (flat > 70) & np.isfinite(p)
    apnd.append(np.sum(p[k] * AREA[k]) / np.sum(AREA[k])); years.append(y)
years = np.array(years)
for lo in range(Y0, Y1 + 1, 10):
    m = (years >= lo) & (years < lo + 10)
    print(f'{lo}-{lo+9}: NH Sep {np.mean(np.array(nhsep)[m]):.2f}  SH Feb {np.mean(np.array(shfeb)[m]):.2f}  July apnd {np.mean(np.array(apnd)[m]):.3f}  July fal {np.mean(np.array(fal)[m]):.3f}')

SURF, TXT1, TXT2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
fig, axes = plt.subplots(2, 1, figsize=(8.5, 6.5), dpi=150, sharex=True); fig.patch.set_facecolor(SURF)
ax = axes[0]; ax.set_facecolor(SURF)
ax.plot(years, nhsep, color='#2a78d6', lw=1.6, label='NH September extent')
ax.plot(years, shfeb, color='#eb6834', lw=1.6, label='SH February extent')
ax.axhline(7.0, color='#2a78d6', ls=':', lw=1); ax.axhline(3.0, color='#eb6834', ls=':', lw=1)
ax.set_ylabel('M km$^2$', color=TXT2, fontsize=9); ax.legend(fontsize=8, frameon=False, loc='center left')
ax.set_title('(a) summer sea-ice extent (dotted: satellite era)', loc='left', fontsize=10, color=TXT1)
ax2 = axes[1]; ax2.set_facecolor(SURF)
ax2.plot(years, apnd, color='#1baf7a', lw=1.6, label='FESOM July pond fraction')
ax2.set_ylabel('pond fraction', color=TXT2, fontsize=9)
ax3 = ax2.twinx(); ax3.plot(years, fal, color='#6b4fbb', lw=1.6, label='OpenIFS July surface albedo')
ax3.set_ylabel('albedo', color=TXT2, fontsize=9)
ax2.set_title('(b) Arctic pack ice in July: FESOM ponds and the albedo OpenIFS uses', loc='left', fontsize=10, color=TXT1)
h1, l1 = ax2.get_legend_handles_labels(); h2, l2 = ax3.get_legend_handles_labels()
ax2.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc='center left')
for a_ in (ax, ax2, ax3):
    a_.axvline(1399.5, color='k', lw=0.8, ls='--')
    for sp in ('top',): a_.spines[sp].set_visible(False)
    a_.tick_params(colors=TXT2, labelsize=8.5)
for a_ in (ax, ax2): a_.grid(color=GRID, lw=0.8)
ax2.set_xlabel('model year (constant 1990 forcing; stronger ponds from 1400)', color=TXT2, fontsize=9)
fig.tight_layout()
out = os.path.join(REPO, 'plots', 'pond_extension_16e1990.png'); fig.savefig(out, facecolor=SURF); print('saved', out)
