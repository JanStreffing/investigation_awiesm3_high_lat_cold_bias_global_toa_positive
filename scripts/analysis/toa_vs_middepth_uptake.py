"""Does lower mid-depth (700-2000 m) ocean heat uptake go with a lower TOA imbalance in decade 4?

From data/coupled_annual_diag.nc.  Decade 4 = years 31-40 of each run (y0+30 .. y0+39).
  toa4     mean net TOA over decade 4 [W/m2]
  mid4     700-2000 m heating rate over decade 4: [H(y0+39) - H(y0+29)] / 10 yr  [W/m2 of Earth]
  mid40    700-2000 m heating rate averaged over years 1-40: [H(y0+39) - H(y0)] / 39 yr
  ohu4     total (all depths) heating rate over decade 4, for the energy-conservation check
Excludes the 1990-forced runs (11P, 11R, 11V, *_1990), the SP-leak runs (11Y, 15A-C), the CMIP7
spin-ups and runs shorter than 40 years.  Prints the table and Pearson r (n, p) for each pair
and writes report/plots/toa_vs_middepth_uptake.png.

Usage:  python3 scripts/analysis/toa_vs_middepth_uptake.py
"""
import os, re
import numpy as np, xarray as xr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AE, SY = 5.101e14, 365.25 * 86400
EXCL = {'11P', '11R', '11V', '11Y', '15A', '15B', '15C'}
ds = xr.open_dataset(os.path.join(REPO, 'data', 'coupled_annual_diag.nc'))
yrs = ds['year'].values.astype(int)

def short(run):
    s = run.replace('Tuning_test_', ''); m = re.match(r'^(\d{2,3}[A-Z]?\d?)', s)
    return m.group(1) if m else s

rows = []
for run in ds['run'].values.astype(str):
    lab = short(run)
    if lab in EXCL or run.endswith('_1990') or run.startswith(('AWIESM7', 'Test_')):
        continue
    t = ds['toa'].sel(run=run).values
    H = {v: ds[v].sel(run=run).values for v in ('ohc_0_100', 'ohc_100_700', 'ohc_700_2000', 'ohc_gt2000')}
    ok = np.isfinite(t)
    if ok.sum() < 40:
        continue
    y0 = yrs[ok][0]; ix = {y: j for j, y in enumerate(yrs)}
    dec4 = [ix[y] for y in range(y0 + 30, y0 + 40) if y in ix]
    need = [y0, y0 + 29, y0 + 39]
    if len(dec4) < 10 or any(y not in ix or not np.isfinite(H['ohc_700_2000'][ix[y]]) for y in need):
        continue
    toa4 = float(np.nanmean(t[dec4]))
    mid = H['ohc_700_2000']
    mid4 = (mid[ix[y0 + 39]] - mid[ix[y0 + 29]]) / (10 * SY * AE)
    mid40 = (mid[ix[y0 + 39]] - mid[ix[y0]]) / (39 * SY * AE)
    tot = sum(H.values())
    ohu4 = (tot[ix[y0 + 39]] - tot[ix[y0 + 29]]) / (10 * SY * AE)
    rows.append((lab, toa4, mid4, mid40, ohu4))

rows.sort(key=lambda r: r[1])
print(f'{"run":<7}{"toa4":>8}{"mid4":>8}{"mid40":>8}{"ohu4":>8}{"mid4/ohu4":>11}')
for lab, a, b, c, d in rows:
    print(f'{lab:<7}{a:8.2f}{b:8.2f}{c:8.2f}{d:8.2f}{(b/d if d else np.nan):11.2f}')
A = np.array([r[1:] for r in rows]); n = len(A)

def pearson(x, y):
    r = np.corrcoef(x, y)[0, 1]
    try:
        from scipy import stats; p = stats.pearsonr(x, y)[1]
    except Exception:
        p = np.nan
    return r, p
print(f'\nn = {n} runs (1850 forcing, >= 40 yr, no SP leak)')
for j, nm in ((1, 'mid4  (700-2000 m rate, decade 4)'), (2, 'mid40 (700-2000 m rate, years 1-40)'),
              (3, 'ohu4  (total ocean uptake, decade 4)')):
    r, p = pearson(A[:, j], A[:, 0]); s = np.polyfit(A[:, j], A[:, 0], 1)[0]
    print(f'  toa4 vs {nm}: r = {r:+.2f}  p = {p:.3f}  slope = {s:+.2f} W/m2 per W/m2')

SURF, TXT1, TXT2, GRID, C = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df', '#2a78d6'
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), dpi=150, sharey=True); fig.patch.set_facecolor(SURF)
for ax, j, xl in ((axes[0], 1, '700-2000 m heating rate, decade 4  [W m$^{-2}$ of Earth]'),
                  (axes[1], 2, '700-2000 m heating rate, years 1-40  [W m$^{-2}$ of Earth]')):
    ax.set_facecolor(SURF)
    ax.scatter(A[:, j], A[:, 0], s=36, color=C, edgecolor=SURF, linewidth=1.2, zorder=3)
    for (lab, *v) in rows:
        ax.annotate(lab, (v[j], v[0]), xytext=(4, 3), textcoords='offset points', fontsize=7.5, color=TXT2)
    s, b = np.polyfit(A[:, j], A[:, 0], 1); xx = np.linspace(A[:, j].min(), A[:, j].max(), 2)
    ax.plot(xx, s * xx + b, color=TXT2, lw=1, ls='--', zorder=2)
    r, p = pearson(A[:, j], A[:, 0])
    ax.text(0.03, 0.95, f'r = {r:+.2f}, p = {p:.3f}, n = {n}', transform=ax.transAxes, fontsize=8.5, color=TXT1, va='top')
    ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=TXT2, labelsize=8.5); ax.set_xlabel(xl, color=TXT2, fontsize=9)
axes[0].set_ylabel('net TOA, decade 4  [W m$^{-2}$]', color=TXT2, fontsize=9)
fig.suptitle('TOA imbalance in decade 4 against mid-depth ocean heat uptake (1850 runs, >= 40 yr)',
             x=0.01, ha='left', fontsize=11, color=TXT1)
fig.tight_layout()
out = os.path.join(REPO, 'report', 'plots', 'toa_vs_middepth_uptake.png'); fig.savefig(out, facecolor=SURF)
print('saved', out)
