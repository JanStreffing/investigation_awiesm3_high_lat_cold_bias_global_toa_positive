"""Ocean heating rate by depth range in 5-yr chunks, every coupled run, as heatmaps.

Reads data/coupled_annual_diag.nc (scripts/analysis/coupled_annual_store.py): annual heat
content by depth range.  Rate for chunk k = [H(y0+5k+5) - H(y0+5k)] / 5 yr,
in W/m2 of Earth area, so it compares directly with net TOA.  Two panels share the rows:
700-2000 m and >2000 m.  Diverging scale centred on 0 (blue = cooling, red = warming).
Also prints the table for all four depth ranges.

Usage:  python3 scripts/figures/ohc_rate_heatmap_all_runs.py
"""
import os, csv, re
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AE, SY, CH, ND = 5.101e14, 365.25 * 86400, 5, 9
FORCED_1990 = {'11P', '11R', '11V'}
SP_LEAK = {'11Y', '15A', '15B', '15C'}
CURRENT = {'15F', '16A', '16B', '16C', '16D'}
COLS = ['J_0_100', 'J_100_700', 'J_700_2000', 'J_gt2000']
NAMES = ['0-100 m', '100-700 m', '700-2000 m', '>2000 m']


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
    dig = re.match(r'^\d+', s)
    dig = dig.group(0) if dig else ''
    num = (int(dig) / 10 if len(dig) == 3 else int(dig)) if dig else 999
    return (1, num, s)


import xarray as xr
DS = xr.open_dataset(os.path.join(REPO, 'data', 'coupled_annual_diag.nc'))
OHC = ['ohc_0_100', 'ohc_100_700', 'ohc_700_2000', 'ohc_gt2000']
H = {}
for run in DS['run'].values.astype(str):
    a = np.stack([DS[v].sel(run=run).values for v in OHC], axis=1)   # (year, 4)
    ok = np.all(np.isfinite(a), axis=1)
    if ok.sum() >= 2:
        H[run] = {int(y): a[j] for j, y in enumerate(DS['year'].values) if ok[j]}
runs = sorted([r for r in H if len(H[r]) >= 2], key=order_key)
RATE = np.full((len(runs), ND, 4), np.nan)
SHORT = np.zeros((len(runs), ND), bool)   # chunk measured over fewer than CH years
for i, run in enumerate(runs):
    y0 = min(H[run])
    for k in range(ND):
        a, b = y0 + CH * k, y0 + CH * (k + 1)
        if a in H[run] and b in H[run]:
            RATE[i, k] = (H[run][b] - H[run][a]) / (CH * SY * AE)
        elif a in H[run]:
            inner = [y for y in H[run] if a < y < b]
            if inner:
                e = max(inner); RATE[i, k] = (H[run][e] - H[run][a]) / ((e - a) * SY * AE); SHORT[i, k] = True

lab = lambda k: f'yr {CH*k+1}-{CH*k+CH}'
for j, nm in enumerate(NAMES):
    print(f'\nocean heating rate {nm} [W/m2 of Earth], {CH}-yr chunks since each run\'s start')
    print(f'  {"run":<18}' + ''.join(f'{lab(k):>10}' for k in range(ND)))
    for i, run in enumerate(runs):
        print(f'  {short(run):<18}' + ''.join((f'{RATE[i,k,j]:9.2f}' + ('*' if SHORT[i, k] else ' ')) if np.isfinite(RATE[i, k, j]) else f'{"":>10}' for k in range(ND)))

BLUE = ['#184f95', '#256abf', '#3987e5', '#6da7ec', '#9ec5f4', '#cde2fb']
MID = '#f0efec'
RED = ['#fbd3d2', '#f5a3a1', '#ec7472', '#e34948', '#c0302f', '#8f1f1f']
cmap = ListedColormap(BLUE + [MID] + RED)
LIM = 1.625
norm = BoundaryNorm(np.linspace(-LIM, LIM, len(BLUE) + len(RED) + 2), cmap.N)
SURF, TXT1, TXT2 = '#fcfcfb', '#0b0b0b', '#52514e'

fig, axes = plt.subplots(1, 2, figsize=(17, 0.26 * len(runs) + 1.8), dpi=150, sharey=True)
fig.patch.set_facecolor(SURF)
for ax, j in zip(axes, (2, 3)):
    ax.set_facecolor(SURF)
    M = RATE[:, :, j]
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, norm=norm, aspect='auto')
    for i in range(len(runs)):
        for k in range(ND):
            if np.isfinite(M[i, k]):
                v = M[i, k]
                ax.text(k, i, f'{v:+.2f}' + ('*' if SHORT[i, k] else ''), ha='center', va='center', fontsize=6.3,
                        color='white' if abs(v) > 0.95 else TXT1)
    ax.set_xticks(range(ND)); ax.set_xticklabels([lab(k) for k in range(ND)], fontsize=7.5, color=TXT2)
    ax.xaxis.tick_top()
    ax.set_xticks(np.arange(-.5, ND, 1), minor=True); ax.set_yticks(np.arange(-.5, len(runs), 1), minor=True)
    ax.grid(which='minor', color=SURF, lw=1.5); ax.tick_params(which='both', length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(f'{NAMES[j]}', loc='left', fontsize=10, color=TXT1, pad=20)
axes[0].set_yticks(range(len(runs)))
labs = axes[0].set_yticklabels([short(r) for r in runs], fontsize=7.5, color=TXT2)
for t, r in zip(labs, runs):
    if short(r).split(' ')[0] in CURRENT:
        t.set_fontweight('bold'); t.set_color(TXT1)
cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.01, ticks=[-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
cb.set_label('ocean heating rate [W m$^{-2}$ of Earth], 5-yr mean', fontsize=8, color=TXT2)
cb.ax.tick_params(labelsize=7, colors=TXT2); cb.outline.set_visible(False)
fig.suptitle('Deep-ocean heat uptake in 5-yr chunks, every coupled run (red = the layer is warming; * = chunk shorter than 5 yr)',
             x=0.01, ha='left', fontsize=11, color=TXT1)
out = os.path.join(REPO, 'report', 'plots', 'ohc_rate_heatmap_all_runs.png')
fig.savefig(out, facecolor=SURF, bbox_inches='tight')
print('\nsaved', out)
