"""Labrador March sea ice across the run history, coloured by the GM/Redi setting.

A: every run decade, March ice in the convection box against global-mean SST.
B: the start-up shock, March ice (5-yr running mean) against years since the run start.
C: the one-change pair PI200_gmR1500 / gmR2500 (same branch, same years).
D: the long runs, 11-yr running mean against years since the start of their line.
Data: data/clim/labrador_ice_vs_gsst_decades.csv, data/clim/lab_hist, data/clim/lab_ice.
GM classes from each run's namelist.oce (K_GM_max / Redi_Kmax); CORE2 AWI-CM3 uses a
different GM formulation (K_GM_max 3000 ramped off below 30-40 km, Ferreira scaling) and HR
is unchecked, so both are grey.
"""
import glob, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
C1000, C1500, C2500, CGREY = '#86b6ef', '#2a78d6', '#0d366b', '#a3a29b'
INK, INK2 = '#1f1f1e', '#6b6a63'
GM2500 = ('16E', 'PI200', 'PICAL', 'LR01', 'LR02', 'LR03')
def gm(e):
    if e.startswith(('CM33', 'CMdev', 'HR')): return 'other'
    if e.startswith(('16C', '16D')): return '1500'
    if e == 'PI200_gmR1500': return '1500'
    if e.startswith(GM2500): return '2500'
    return '1000'
COL = {'1000': C1000, '1500': C1500, '2500': C2500, 'other': CGREY}
LAB = {'1000': 'K_GM_max 1000, Redi 0', '1500': 'K_GM_max 1500, Redi 1000',
       '2500': 'K_GM_max 2500, Redi 1000', 'other': 'AWI-CM3 v3.3 / HR\n(other GM formulation)'}
# annual March box ice
h = pd.concat(pd.read_csv(f, names=['exp', 'year', 'month', 'region', 'var', 'v']) for f in glob.glob('data/clim/lab_hist/*.csv'))
h = h[(h.month == 3) & (h.region == 'interior') & (h['var'] == 'a_ice')][['exp', 'year', 'v']]
l = pd.concat(pd.read_csv(f, names=['exp', 'year', 'month', 'region', 'v']) for f in glob.glob('data/clim/lab_ice/lab_ice_*.csv'))
l = l[(l.month == 3) & (l.region == 'interior')][['exp', 'year', 'v']]
A = pd.concat([h[~h.exp.isin(l.exp.unique())], l]).drop_duplicates(['exp', 'year'], keep='last')
def s(e, a=None, b=None):
    x = A[A.exp == e].set_index('year').v.sort_index()
    return x.loc[a:b] if a is not None else x
D = pd.read_csv('data/clim/labrador_ice_vs_gsst_decades.csv')
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#c9c8c0', 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.titlesize': 10, 'axes.titleweight': 'bold'})
fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.5))
for a in ax.flat:
    a.grid(color='#e6e5df', lw=0.6); a.set_axisbelow(True)
    for sp in ('top', 'right'): a.spines[sp].set_visible(False)
OBS = 0.16
# A
a = ax[0, 0]
# first two decades of a run that starts from a new ocean state = start-up shock (hollow)
first = D.groupby('exp').y0.transform('min'); early = (D.y0 - first) < 20
for k in ('other', '1000', '1500', '2500'):
    m = D.exp.map(gm) == k
    d = D[m & ~early]; e = D[m & early]
    a.scatter(d.gsst, d.ice_box, s=26 if k != 'other' else 18, color=COL[k], edgecolor='white', lw=0.6,
              label=f'{LAB[k]}  ({D[m].exp.nunique()} runs)', zorder=3 if k == '2500' else 2)
    a.scatter(e.gsst, e.ice_box, s=26 if k != 'other' else 18, facecolor='none', edgecolor=COL[k], lw=1.1,
              zorder=3 if k == '2500' else 2)
a.scatter([], [], s=26, facecolor='none', edgecolor=INK2, lw=1.1, label='hollow: first 20 yr of a run (start-up)')
a.axhline(OBS, color=INK2, lw=1, ls='--'); a.text(20.45, OBS + 0.012, 'HadISST2 1979-2008', ha='right', fontsize=8, color=INK2)
a.set_xlabel('global-mean SST, decade mean (°C)'); a.set_ylabel('March sea-ice fraction, convection box')
a.set_title('A  Every run decade: beyond start-up, high ice only with GM 2500', loc='left')
a.legend(fontsize=7.5, frameon=False, loc='upper right', bbox_to_anchor=(1.0, 0.98))
# B: start-up shock
a = ax[0, 1]
runs = [('11X', '1000'), ('15F', '1000'), ('TT06_Baseline', '1000'), ('16C', '1500'), ('16D', '1500'),
        ('16E', '2500'), ('LR01_CMIP7hist', '2500'), ('PICAL', '2500')]
for e, k in runs:
    x = s(e); 
    if e == 'PICAL':   # PICAL lineage start, first 60 years
        x = x.loc[1850:1909]
    x = x.iloc[:60]; t = np.arange(len(x))
    y = x.rolling(5, center=True, min_periods=3).mean()
    a.plot(t, y.values, color=COL[k], lw=2)
    name = {'TT06_Baseline': '06 Baseline', 'LR01_CMIP7hist': 'LR01 (hist.)'}.get(e, e)
    dy = {'11X': 0.035, '16D': 0.0, '16C': -0.035, '15F': 0.02, 'TT06_Baseline': -0.01}.get(e, 0)
    a.text(t[-1] + 1, y.values[-1] + dy, name, fontsize=7.5, color=INK, va='center')
a.axhline(OBS, color=INK2, lw=1, ls='--')
a.set_xlim(0, 70); a.set_xlabel('years since run start'); a.set_ylabel('March sea-ice fraction, 5-yr running mean')
a.set_title('B  Start-up shock grows with GM: 1000 < 1500 < 2500', loc='left')
# C: clean pair
a = ax[1, 0]
p = s('PI200', 1555, 1580)
a.plot(p.index, p.values, color=CGREY, lw=1.6, marker='o', ms=3.5, label='PI200 (parent, GM 2500, resolution scaling)')
for e, k, lab in (('PI200_gmR1500', '1500', 'gmR1500: Rossby scaling, K_GM_max 1500'),
                  ('PI200_gmR2500', '2500', 'gmR2500: Rossby scaling, K_GM_max 2500')):
    x = s(e); a.plot(x.index, x.values, color=COL[k], lw=2, marker='o', ms=4, label=lab)
    a.hlines(x.mean(), x.index.min(), x.index.max(), color=COL[k], lw=1, ls=':')
    a.text(x.index.max() + 0.6, x.mean(), f'mean {x.mean():.2f}', fontsize=8, color=INK, va='center')
a.axhline(OBS, color=INK2, lw=1, ls='--'); a.axvline(1579.5, color='#c9c8c0', lw=1)
a.text(1579.8, 0.97, 'branch 1580', fontsize=7.5, color=INK2, va='top')
a.set_xlim(1554, 1604); a.set_ylim(0, 1)
a.set_xlabel('model year'); a.set_ylabel('March sea-ice fraction, convection box')
a.set_title('C  One change, same years: GM 2500 0.68, GM 1500 0.28', loc='left')
a.legend(fontsize=7.5, frameon=False, loc='upper left')
# D: long runs
a = ax[1, 1]
long = [('AWIESM720_SPINUP', '1000', 1350, 'AWIESM700/720 spin-ups'), ('AWIESM700_SPINUP', '1000', 1350, None),
        ('16E_1990', '2500', 1350, None), ('PI200', '2500', 1350, '16E_1990 -> PI200'),
        ('PICAL', '2500', 1850, None), ('PICAL_momixoff', '2500', 1850, None), ('PICAL_ccnice', '2500', 1850, 'PICAL -> ccnice'),
        ('PICAL_momixoff_2040', '2500', 1850, 'momixoff continued')]
for e, k, t0, lab in long:
    x = s(e).rolling(11, center=True, min_periods=6).mean()
    ls = '-' if e != 'PICAL_momixoff_2040' else (0, (4, 2))
    a.plot(x.index - t0, x.values, color=COL[k], lw=1.8, ls=ls)
    if lab:
        a.text(x.index[-1] - t0 + 4, x.values[-1], lab, fontsize=7.5, color=INK, va='center')
a.axhline(OBS, color=INK2, lw=1, ls='--')
a.set_xlim(0, 590); a.set_ylim(0, 1)
a.set_xlabel('years since the line started'); a.set_ylabel('March sea-ice fraction, 11-yr running mean')
a.text(300, 0.24, 'caution: these drift to 19.2 °C global SST,\n~1 K warmer than the PICAL line (18.2-18.4)', fontsize=7.5, color=INK2)
a.set_title('D  Long runs: GM 1000 never caps in 500 yr; GM 2500 lines do', loc='left')
fig.suptitle('Labrador Sea March ice across the run history, by GM / Redi setting', fontsize=12, fontweight='bold', color=INK, x=0.01, ha='left')
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig('plots/labrador_ice_gm_history.png', dpi=150)
