"""Fox-Kemper MLE from clean pairs: 15F -> 16B (10 yr) and 16C -> 16D (GM 1500 base, ~25 yr).

(a) arm minus control for heating rate by depth range, total ocean uptake and net TOA
    [W/m2 of Earth], mean over common years, error bar 1.96 x s.e. of the annual difference;
(b) the same for the AMOC indices [Sv];
(c) winter mixed-layer depth bias against WOA18 (deepest month, WOA criterion), control -> arm
    per band.  Values from scripts/analysis/mld_baseline_vs_woa18.py: 15F/16B 1355-59,
    16C/16D 1365-69 (tables of 2026-09-10/11).
Usage:  python3 scripts/figures/mle_pair_effects.py            (both pairs)
        PAIRS=15F:16B python3 scripts/figures/mle_pair_effects.py   (one pair)
"""
import os
import numpy as np, xarray as xr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AE, SY = 5.101e14, 365.25 * 86400
C = xr.open_dataset(os.path.join(REPO, 'data', 'coupled_annual_diag.nc'))
A = xr.open_dataset(os.path.join(REPO, 'data', 'amoc_annual_diag.nc'))
ALL = {'15F:16B': ('15F', '16B', '15F → 16B (use_mle only, 10 yr)', '#2a78d6'),
       '16C:16D': ('16C', '16D', '16C → 16D (use_mle only, 25 yr)', '#eb6834')}
SEL = os.environ.get('PAIRS', '15F:16B,16C:16D').split(',')
PAIRS = [ALL[k] for k in SEL]
DEP = [('0-100 m', 'ohc_0_100'), ('100-700 m', 'ohc_100_700'), ('700-2000 m', 'ohc_700_2000'), ('>2000 m', 'ohc_gt2000')]

def diffs(ctl, arm):
    out = {}
    y = C['year'].values.astype(int)
    tot = {r: 0 for r in (ctl, arm)}
    for lab, v in DEP + [('total uptake', None), ('net TOA', 'toa')]:
        if v is None:
            a = sum(np.diff(C[vv].sel(run=arm).values) for _, vv in DEP) / (SY * AE)
            c = sum(np.diff(C[vv].sel(run=ctl).values) for _, vv in DEP) / (SY * AE)
        elif v == 'toa':
            a, c = C[v].sel(run=arm).values, C[v].sel(run=ctl).values
        else:
            a, c = np.diff(C[v].sel(run=arm).values) / (SY * AE), np.diff(C[v].sel(run=ctl).values) / (SY * AE)
        d = (a - c)[np.isfinite(a - c)]
        out[lab] = (d.mean(), 1.96 * d.std(ddof=1) / np.sqrt(d.size))
    for lab, v in (('AMOC 26.5N', 'amoc26'), ('AMOC 40-60N', 'amoc4060')):
        d = (A[v].sel(run=arm).values - A[v].sel(run=ctl).values); d = d[np.isfinite(d)]
        out[lab] = (d.mean(), 1.96 * d.std(ddof=1) / np.sqrt(d.size))
    return out

MLD = {  # deepest-month bias vs WOA18 [m], band order north -> south
    'bands': ['60-90N', '45-60N', '30-45N', '30S-30N', '45-30S', '60-45S', '90-60S'],
    '15F': [-16.4, -45.7, 50.4, 6.5, 34.7, -39.5, 593.6], '16B': [19.1, -12.7, 28.3, -3.7, 2.1, -102.9, 606.2],
    '16C': [42.0, -3.2, 41.5, 0.1, 50.4, -23.7, 860.5], '16D': [18.6, -51.2, 17.6, -8.4, 6.1, -73.0, 742.5]}

SURF, TXT1, TXT2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
fig = plt.figure(figsize=(13, 7.6), dpi=150); fig.patch.set_facecolor(SURF)
gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[2.2, 1], hspace=0.45, wspace=0.28)
axa, axb, axc = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])
D = [diffs(c, a) for c, a, _, _ in PAIRS]

def forest(ax, labels, unit, title):
    ax.set_facecolor(SURF); yv = np.arange(len(labels))[::-1]
    for k, ((c, a, lab, col), d) in enumerate(zip(PAIRS, D)):
        off = (0.14 if k == 0 else -0.14) if len(PAIRS) > 1 else 0.0
        m = [d[l][0] for l in labels]; e = [d[l][1] for l in labels]
        ax.errorbar(m, yv + off, xerr=e, fmt='o', ms=6, color=col, ecolor=col, elinewidth=1.6, capsize=3,
                    mec=SURF, mew=1, label=lab, zorder=3)
    ax.axvline(0, color=TXT2, lw=0.9, zorder=1)
    ax.set_yticks(yv); ax.set_yticklabels(labels, fontsize=9, color=TXT1)
    ax.grid(axis='x', color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color(GRID); ax.tick_params(colors=TXT2, labelsize=8.5, length=0)
    ax.set_xlabel(unit, color=TXT2, fontsize=9); ax.set_title(title, loc='left', fontsize=10, color=TXT1)

forest(axa, ['0-100 m', '100-700 m', '700-2000 m', '>2000 m', 'total uptake', 'net TOA'],
       'MLE arm minus control  [W m$^{-2}$ of Earth]', '(a) Ocean heating rate by depth, and net TOA')
h, l = axa.get_legend_handles_labels()
fig.legend(h, l, loc='upper left', bbox_to_anchor=(0.01, 0.955), ncol=2, fontsize=9, frameon=False)
forest(axb, ['AMOC 26.5N', 'AMOC 40-60N'], 'MLE arm minus control  [Sv]', '(b) AMOC')

axc.set_facecolor(SURF); yb = np.arange(len(MLD['bands']))[::-1]
for k, (c, a, lab, col) in enumerate(PAIRS):
    off = (0.16 if k == 0 else -0.16) if len(PAIRS) > 1 else 0.0
    for i, b in enumerate(MLD['bands']):
        x0, x1 = MLD[c][i], MLD[a][i]
        axc.annotate('', xy=(x1, yb[i] + off), xytext=(x0, yb[i] + off),
                     arrowprops=dict(arrowstyle='-|>', color=col, lw=1.6, mutation_scale=10), zorder=3)
        axc.plot([x0], [yb[i] + off], 'o', ms=4.5, color=SURF, mec=col, mew=1.4, zorder=4)
axc.axvline(0, color=TXT2, lw=0.9)
axc.set_xscale('symlog', linthresh=100); axc.set_xlim(-160, 1200)
axc.set_xticks([-100, -50, 0, 50, 100, 500, 1000]); axc.set_xticklabels(['-100', '-50', '0', '50', '100', '500', '1000'])
axc.set_yticks(yb); axc.set_yticklabels(MLD['bands'], fontsize=9, color=TXT1)
axc.grid(axis='x', color=GRID, lw=0.8); axc.set_axisbelow(True)
for s in ('top', 'right', 'left'): axc.spines[s].set_visible(False)
axc.spines['bottom'].set_color(GRID); axc.tick_params(colors=TXT2, labelsize=8.5, length=0)
axc.set_xlabel('winter mixed-layer depth bias vs WOA18 [m]\n(open dot = control, arrow head = MLE arm; symlog)', color=TXT2, fontsize=9)
axc.set_title('(c) Winter mixed layer, control → MLE', loc='left', fontsize=10, color=TXT1)
fig.suptitle('Fox–Kemper MLE: ' + ('clean pairs' if len(PAIRS) > 1 else PAIRS[0][2]) + ', arm minus control (error bars 95 %)', x=0.01, ha='left', fontsize=12, color=TXT1)
out = os.path.join(REPO, 'report', 'plots', 'mle_pair_effects.png' if len(PAIRS) > 1 else f'mle_pair_effects_{PAIRS[0][0]}_{PAIRS[0][1]}.png'); fig.savefig(out, facecolor=SURF, bbox_inches='tight')
print('saved', out)
for (c, a, lab, _), d in zip(PAIRS, D):
    print(lab); [print(f'  {k:<13} {v[0]:+.3f} +- {v[1]:.3f}') for k, v in d.items()]
