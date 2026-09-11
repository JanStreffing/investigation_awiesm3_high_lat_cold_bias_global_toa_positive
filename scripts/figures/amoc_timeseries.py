"""AMOC time series for 15F and the 16 series against the older long coupled runs.

Reads data/amoc_annual_diag.nc (scripts/analysis/amoc_annual_store.py).  Two panels share
the x axis: the maximum at 26.5N and over 40-60N, both below 500 m.  Solid = Eulerian (w),
dashed = residual (w + bolus_w, adds the GM eddy-induced part).  Older runs in grey.
Also prints decadal means of both indices.

Usage:  python3 scripts/figures/amoc_timeseries.py
"""
import os
import numpy as np, xarray as xr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ds = xr.open_dataset(os.path.join(REPO, 'data', 'amoc_annual_diag.nc'))
yrs = ds['year'].values.astype(int)
MAIN = [('15F', '#2a78d6'), ('16A', '#eb6834'), ('16B', '#1baf7a'), ('16C', '#eda100'), ('16D', '#e87ba4'),
        ('16E', '#6b4fbb')]
CTX = [('11X', '-'), ('11W', '--'), ('11E', ':'), ('11G', '-.'), ('11I', (0, (1, 3)))]
runs = list(ds['run'].values.astype(str))

print('decadal mean AMOC [Sv]  (26.5N Eulerian / residual | 40-60N Eulerian / residual)')
for r in [x for x, _ in CTX] + [x for x, _ in MAIN]:
    if r not in runs:
        continue
    s = ds.sel(run=r); cells = []
    for lo in range(1350, 1400, 10):
        k = (yrs >= lo) & (yrs < lo + 10)
        v = [np.nanmean(s[n].values[k]) if np.isfinite(s[n].values[k]).any() else np.nan
             for n in ('amoc26', 'amoc26_res', 'amoc4060', 'amoc4060_res')]
        cells.append('      --      ' if np.isnan(v[0]) else f'{v[0]:4.1f}/{v[1]:4.1f}|{v[2]:4.1f}/{v[3]:4.1f}')
    print(f'  {r:<5}' + '  '.join(cells))

SURF, TXT1, TXT2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), dpi=150, sharex=True); fig.patch.set_facecolor(SURF)
for ax, (e, rsd, title) in zip(axes, (('amoc26', 'amoc26_res', 'AMOC at 26.5N (max below 500 m)'),
                                      ('amoc4060', 'amoc4060_res', 'AMOC maximum over 40-60N, below 500 m'))):
    ax.set_facecolor(SURF)
    for r, ls in CTX:
        if r not in runs: continue
        v = ds[e].sel(run=r).values; k = np.isfinite(v)
        ax.plot(yrs[k], v[k], color='#8a8985', lw=1.0, ls=ls, label=r, zorder=2)
        ax.annotate(r, (yrs[k][-1], v[k][-1]), xytext=(4, 0), textcoords='offset points', fontsize=7.5, color=TXT2, va='center')
    for r, c in MAIN:
        if r not in runs: continue
        v = ds[e].sel(run=r).values; k = np.isfinite(v)
        ax.plot(yrs[k], v[k], color=c, lw=2.0 if r == '15F' else 1.5, label=r, zorder=4 if r == '15F' else 3)
        vr = ds[rsd].sel(run=r).values; kr = np.isfinite(vr)
        ax.plot(yrs[kr], vr[kr], color=c, lw=1.0, ls='--', alpha=0.8, zorder=3)
        ax.annotate(r, (yrs[k][-1], v[k][-1]), xytext=(4, 0), textcoords='offset points', fontsize=8.5,
                    color=TXT1, fontweight='bold' if r == '15F' else 'normal', va='center')
    ax.set_title(title, loc='left', fontsize=10, color=TXT1)
    ax.grid(axis='y', color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=TXT2, labelsize=8.5); ax.set_ylabel('Sv', color=TXT2, fontsize=9)
axes[1].set_xlabel('model year', color=TXT2, fontsize=9)
axes[0].legend(loc='upper right', fontsize=8, frameon=False, ncol=5)
fig.suptitle('AMOC, 15F and the 16 series against the older long runs (solid Eulerian, dashed residual incl. GM)',
             x=0.01, ha='left', fontsize=11, color=TXT1)
fig.tight_layout()
out = os.path.join(REPO, 'report', 'plots', 'amoc_timeseries.png'); fig.savefig(out, facecolor=SURF)
print('saved', out)
