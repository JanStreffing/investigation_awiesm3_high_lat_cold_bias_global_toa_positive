"""Decadal-mean global net TOA for every coupled run, as a heatmap (runs x decades since start).

Reads data/coupled_annual_diag.nc (scripts/analysis/coupled_annual_store.py), runs with >= 10 years.  Columns are decades
counted from each run's first year, so the 1850-dated piControl lines up with the 1350-dated
tuning arms.  Diverging scale centred on 0: blue = losing energy, red = gaining.
Rows marked (1990) ran with constant 1990 forcing, (SP leak) on the unfixed single-precision ocean, (piCtl) is the concentration-driven
piControl; the current arms (15F, 16A-C) are in bold.  Also prints the table.

Usage:  python3 scripts/figures/toa_decadal_heatmap_all_runs.py            (10-yr chunks)
        CHUNK=5 python3 scripts/figures/toa_decadal_heatmap_all_runs.py    (5-yr chunks)
"""
import os, csv, re
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import xarray as xr
DS = xr.open_dataset(os.path.join(REPO, 'data', 'coupled_annual_diag.nc'))
FORCED_1990 = {'11P', '11R', '11V'}
SP_LEAK = {'11Y', '15A', '15B', '15C'}       # single-precision FESOM before the FCT fix
EXCLUDE = ('AWIESM7',)   # CMIP7 spin-ups: TOA comes out at ~+10 W/m2, a different output convention
CURRENT = {'15F', '16A', '16B', '16C', '16D'}


def short(run):
    s = run.replace('Tuning_test_', '')
    m = re.match(r'^(\d{2,3}[A-Z]?\d?)', s)
    lab = m.group(1) if m else s
    if run.startswith('Test_piControl'):
        lab = 'piCtl ' + s.replace('Test_piControl_TCO095_CORE3_', '')[:22]
    if lab in FORCED_1990:
        lab += ' (1990)'
    if lab in SP_LEAK:
        lab += ' (SP leak)'
    return lab


def order_key(run):
    s = short(run)
    if s.startswith('piCtl'):
        return (2, s)
    m = re.match(r'^(\d{2,3})([A-Z]?)(\d?)', s)
    dig = re.match(r'^\d+', s).group(0) if m else ''
    num = (int(dig) / 10 if len(dig) == 3 else int(dig)) if dig else 999   # 070/080/110 sit with 07/08/11
    return (1, num, s)


by = {}
for run in DS['run'].values.astype(str):
    if run.startswith(EXCLUDE):
        continue
    v = DS['toa'].sel(run=run).values; ok = np.isfinite(v)
    if ok.sum() >= 10:
        by[run] = list(zip(DS['year'].values[ok].astype(int), v[ok]))
runs = sorted(by, key=order_key)
CH = int(os.environ.get('CHUNK', 10))   # chunk length in years
ND = 50 // CH
M = np.full((len(runs), ND), np.nan)
N = np.zeros((len(runs), ND), int)
for i, run in enumerate(runs):
    yv = sorted(by[run]); y0 = yv[0][0]
    for y, v in yv:
        k = (y - y0) // CH
        if k < ND:
            M[i, k] = v if np.isnan(M[i, k]) else M[i, k] + v; N[i, k] += 1
M = np.where(N > 0, M / np.maximum(N, 1), np.nan)

print(f'{CH}-yr mean global net TOA [W/m2], chunks since each run\'s start')
print(f'  {"run":<26}' + ''.join(f'{"yr " + str(CH*k+1) + "-" + str(CH*k+CH):>11}' for k in range(ND)))
for i, run in enumerate(runs):
    cells = ''.join(f'{M[i,k]:7.2f} ({N[i,k]:>2})' if N[i, k] else f'{"":>11}' for k in range(ND))
    print(f'  {short(run):<26}' + cells)

# diverging scale: blue arm = reference palette blue ramp, red arm built around its red pole
BLUE = ['#184f95', '#256abf', '#3987e5', '#6da7ec', '#9ec5f4', '#cde2fb']
MID = '#f0efec'
RED = ['#fbd3d2', '#f5a3a1', '#ec7472', '#e34948', '#c0302f', '#8f1f1f']
cmap = ListedColormap(BLUE + [MID] + RED)
LIM = 1.625
bounds = np.linspace(-LIM, LIM, len(BLUE) + len(RED) + 2)
norm = BoundaryNorm(bounds, cmap.N)

SURF, TXT1, TXT2 = '#fcfcfb', '#0b0b0b', '#52514e'
fig, ax = plt.subplots(figsize=(8.2 if ND <= 5 else 13.5, 0.26 * len(runs) + 1.6), dpi=150)
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, norm=norm, aspect='auto')
for i in range(len(runs)):
    for k in range(ND):
        if N[i, k]:
            v = M[i, k]; dark = abs(v) > 0.95
            ax.text(k, i, f'{v:+.2f}' + ('' if N[i, k] == CH else f' ({N[i,k]})'), ha='center',
                    va='center', fontsize=7 if ND <= 5 else 6.5, color='white' if dark else TXT1)
ax.set_xticks(range(ND)); ax.set_xticklabels([f'yr {CH*k+1}-{CH*k+CH}' for k in range(ND)], fontsize=8, color=TXT2)
ax.xaxis.tick_top()
ax.set_yticks(range(len(runs)))
labs = ax.set_yticklabels([short(r) for r in runs], fontsize=7.5, color=TXT2)
for t, r in zip(labs, runs):
    if short(r).split(' ')[0] in CURRENT:
        t.set_fontweight('bold'); t.set_color(TXT1)
ax.set_xticks(np.arange(-.5, ND, 1), minor=True); ax.set_yticks(np.arange(-.5, len(runs), 1), minor=True)
ax.grid(which='minor', color=SURF, lw=1.5); ax.tick_params(which='minor', length=0)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, ticks=[-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
cb.set_label(f'global net TOA [W m$^{{-2}}$], {CH}-yr mean', fontsize=8, color=TXT2)
cb.ax.tick_params(labelsize=7, colors=TXT2); cb.outline.set_visible(False)
ax.set_title(f'Global net TOA in {CH}-yr chunks, every coupled run (blue = losing energy, red = gaining)',
             loc='left', fontsize=10, color=TXT1, pad=22)
fig.tight_layout()
out = os.path.join(REPO, 'report', 'plots', 'toa_decadal_heatmap_all_runs.png' if CH == 10
                   else f'toa_{CH}yr_heatmap_all_runs.png')
fig.savefig(out, facecolor=SURF)
print('saved', out)
