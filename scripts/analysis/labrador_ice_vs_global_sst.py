"""Labrador March sea ice against global temperature, across every run of the campaign and
before, so that runs that are simply too cold or too warm do not decide the question.

Inputs: data/clim/lab_hist/*.csv (labrador_any_mesh.py, levante runs),
data/clim/lab_ice/{lab_ice,mld}_*.csv (PICAL lineage + albedo arms, CORE3 exact areas),
data/clim/lab_sst/*.csv (global_sst_any_mesh.py). Unit: decade means within each run
(first decade of a run dropped if the run is a branch with a cold-started atmosphere is NOT
attempted; all decades kept, flagged by start year). Output: table + figure.
"""
import glob, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
rows = []
for f in glob.glob('data/clim/lab_hist/*.csv'):
    d = pd.read_csv(f, names=['exp', 'year', 'month', 'region', 'var', 'v']); rows.append(d)
lab = pd.concat(rows) if rows else pd.DataFrame()
# PICAL lineage and albedo arms (exact CORE3 areas) override the levante box extraction
for pat, var in (('data/clim/lab_ice/lab_ice_*.csv', 'a_ice'), ('data/clim/lab_ice/mld_*.csv', 'MLD2')):
    d = pd.concat(pd.read_csv(f, names=['exp', 'year', 'month', 'region', 'v']) for f in glob.glob(pat))
    d['var'] = var
    lab = pd.concat([lab[~lab.exp.isin(d.exp.unique()) | (lab['var'] != var)], d])
lab = lab.drop_duplicates(['exp', 'year', 'month', 'region', 'var'], keep='last')
sst = pd.concat(pd.read_csv(f, names=['exp', 'year', 'gsst', 'nh50', 'sh50'])
                for f in glob.glob('data/clim/lab_sst/*.csv'))
sst = sst.drop_duplicates(['exp', 'year'], keep='last')
mar = lab[lab.month == 3].pivot_table(index=['exp', 'year'], columns=['var', 'region'], values='v')
mar.columns = [f'{a}_{b}' for a, b in mar.columns]
t = mar.join(sst.set_index(['exp', 'year']), how='inner').reset_index()
t['dec'] = t.groupby('exp').year.transform(lambda y: (y - y.min()) // 10)
D = t.groupby(['exp', 'dec']).agg(y0=('year', 'min'), y1=('year', 'max'), n=('year', 'size'),
                                  gsst=('gsst', 'mean'), nh50=('nh50', 'mean'),
                                  ice_box=('a_ice_interior', 'mean'), ice_basin=('a_ice_labsea', 'mean'),
                                  mld_box=('MLD2_interior', 'mean')).reset_index()
D = D[D.n >= 5]
D.to_csv('data/clim/labrador_ice_vs_gsst_decades.csv', index=False)
# Relation across all runs: March box ice vs global SST (robust: binned medians)
ok = D[['gsst', 'ice_box']].notna().all(axis=1)
fit = np.polyfit(D.gsst[ok], D.ice_box[ok], 1)
D['resid'] = D.ice_box - np.polyval(fit, D.gsst)
pd.set_option('display.width', 200); pd.set_option('display.max_rows', 400)
print(f'all decades: {len(D)} from {D.exp.nunique()} runs; linear fit ice_box = {fit[0]:+.3f}*gSST {fit[1]:+.2f}')
print('\nper run (decade means averaged): global SST, NH>50 SST, March ice box/basin, March MLD box, residual vs fit')
R = D.groupby('exp').agg(y0=('y0', 'min'), y1=('y1', 'max'), gsst=('gsst', 'mean'), nh50=('nh50', 'mean'),
                         ice_box=('ice_box', 'mean'), ice_basin=('ice_basin', 'mean'), mld_box=('mld_box', 'mean'),
                         resid=('resid', 'mean')).sort_values('gsst')
print(R.round(3).to_string())
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
fam = lambda e: ('AWI-CM3 / CORE2' if e.startswith(('CM33', 'CMdev', 'TUNE42', 'momix', 'LR_melt')) else
                 'HR' if e.startswith('HR') else
                 'AWIESM7x0 spin-ups' if e.startswith('AWIESM7') else
                 'tuning 06-11V (old CORE3)' if e.startswith(('TT', '11V', 'noLPJG', 'si10y')) else
                 'tuning 11X-16E, PI200' if e.startswith(('11X', '13A', '15F', '16', 'PI200')) else
                 'LR CMIP7 hist' if e.startswith('LR0') else 'PICAL family')
cols = dict(zip(['AWI-CM3 / CORE2', 'HR', 'AWIESM7x0 spin-ups', 'tuning 06-11V (old CORE3)',
                 'tuning 11X-16E, PI200', 'LR CMIP7 hist', 'PICAL family'],
                ['#999999', '#8c564b', '#bcbd22', '#1f77b4', '#ff7f0e', '#17becf', '#d62728']))
for k, (x, ylab) in enumerate((('ice_box', 'March sea-ice fraction, convection box 56-62N 60-50W'),
                               ('ice_basin', 'March sea-ice fraction, Labrador basin 52-66N 65-45W'))):
    for fname, c in cols.items():
        s = D[D.exp.map(fam) == fname]
        ax[k].scatter(s.gsst, s[x], s=14, color=c, alpha=0.7, label=fname)
    ax[k].axhline({'ice_box': 0.16, 'ice_basin': 0.37}[x], color='k', lw=0.8, ls='--', label='HadISST2 1979-2008')
    ax[k].set_xlabel('global-mean SST, decade mean (degC)'); ax[k].set_ylabel(ylab); ax[k].grid(alpha=0.3)
ax[0].legend(fontsize=7)
fig.suptitle('Labrador March sea ice against global temperature, one dot per run decade', fontsize=11)
fig.tight_layout(); fig.savefig('plots/labrador_ice_vs_global_sst.png', dpi=140)
