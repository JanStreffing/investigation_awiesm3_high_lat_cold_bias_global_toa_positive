"""Labrador Sea ice fraction along the tuning lineage, one running year counter.

PICAL 1850-1919 -> PICAL_momixoff 1920-1939 -> PICAL_ccnice 1940-2099 -> PICAL_crunveg
2100-2119 -> PICAL_crunveg_nx1800 2120-.  Each run contributes the years from its own
start to the next branch-off. Input: data/clim/lab_ice/*.csv from
scripts/analysis/labrador_ice_lineage.py (levante for the first three, albedo after).
"""
import glob, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
LIN = [('PICAL', 1850, 1919), ('PICAL_momixoff', 1920, 1939), ('PICAL_ccnice', 1940, 2099),
       ('PICAL_crunveg', 2100, 2119), ('PICAL_crunveg_nx1800', 2120, 2200)]
COL = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3']
df = pd.concat(pd.read_csv(f, names=['exp', 'year', 'month', 'region', 'aice'])
               for f in sorted(glob.glob('data/clim/lab_ice/lab_ice_*.csv')))
df = pd.concat(df[(df.exp == e) & (df.year >= a) & (df.year <= b)] for e, a, b in LIN)
y0 = LIN[0][1]
fig, axs = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True)
for ax, reg, title in zip(axs, ('labsea', 'interior'),
                          ('Labrador Sea basin, 52-66N 65-45W', 'Convection box, 56-62N 60-50W')):
    d = df[df.region == reg]
    ann = d.groupby(['exp', 'year']).aice.mean().reset_index()
    mar = d[d.month == 3].groupby(['exp', 'year']).aice.mean().reset_index()
    for (e, a, b), c in zip(LIN, COL):
        for s, lw, al, lab in ((mar, 0.9, 0.55, None), (ann, 1.8, 1.0, e)):
            x = s[s.exp == e]
            if len(x): ax.plot(x.year - y0 + 1, x.aice, color=c, lw=lw, alpha=al, label=lab)
    for e, a, b in LIN[1:]:
        if a - y0 + 1 <= df.year.max() - y0 + 1:
            ax.axvline(a - y0 + 0.5, color='0.4', lw=0.8, ls=':')
    ax.set_ylabel('sea-ice fraction'); ax.set_title(title, loc='left', fontsize=10)
    ax.set_ylim(bottom=0); ax.grid(alpha=0.3)
axs[0].legend(ncol=5, fontsize=8, loc='upper left', title='thick: annual mean   thin: March', title_fontsize=8)
for e, a, b in LIN[1:]:
    if a <= df.year.max():
        axs[0].annotate(f'{a}', (a - y0 + 0.5, axs[0].get_ylim()[1]), fontsize=7, ha='center', va='bottom', color='0.3')
axs[1].set_xlabel(f'year of the lineage (1 = {y0})')
fig.tight_layout(); fig.savefig('plots/labrador_ice_lineage.png', dpi=150)
print(df.groupby('exp').year.agg(['min', 'max', 'nunique']))
