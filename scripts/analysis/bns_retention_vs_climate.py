"""Is the larch surviving because of the new LPJ-GUESS code, or because the climate is warmer?

WHY THIS EXISTS.  PICAL_crunveg (2100+, merged canopy code) keeps 81 % of its Siberian
BNS after 19 years, where 080a and 11G (same CRUNCEP start, pre-merge binary) keep 34 %
and 20 %.  But crunveg's atmosphere is also warmer (v3.5 tuning and a denser inherited
forest), so the run-mean comparison cannot separate code from climate.  Here each cell is
its own sample: BNS retention (year 19 / year 0) against that cell's own climate over the
same 20 years.  If crunveg retains more at MATCHED climate, the code is doing it; if its
cells simply sit further right on a shared curve, the climate is.

The competition frame (skill bns-not-gate-limited): retention is read against GDD5 and
against the light at the grass layer, not against a gate.  Only cells that start with
real larch (BNS FPC >= 0.05 in year 0) enter.

Usage:  bns_retention_vs_climate.py <out png> <tag>=<root>:<y0> [...]
"""
import glob
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SIB = (55.0, 75.0, 60.0, 180.0)
NYR = 20
COL = {'crunveg': '#eb6834', '080a': '#1baf7a', '11G': '#eda100', 'ccnice': '#2a78d6'}


def load(root, tab, cols):
    frames = []
    for fn in sorted(glob.glob(f'{root}/outdata/lpj_guess/*/run1/{tab}.out')):
        hdr = open(fn).readline().split()
        frames.append(pd.read_csv(fn, sep=r'\s+', usecols=['Lon', 'Lat', 'Year'] + cols))
    d = pd.concat(frames).drop_duplicates(['Lon', 'Lat', 'Year'], keep='last')
    lo = d.Lon % 360
    return d[(d.Lat >= SIB[0]) & (d.Lat <= SIB[1]) & (lo >= SIB[2]) & (lo <= SIB[3])]


def cells(root, y0):
    f = load(root, 'fpc', ['BNS', 'TREEFPC', 'GRASSFPC'])
    e = load(root, 'est_limits', ['agdd5', 'grassPAR', 'mTmin20'])
    dn = load(root, 'dens', ['BNS']).rename(columns={'BNS': 'densBNS'})
    d = f.merge(e, on=['Lon', 'Lat', 'Year']).merge(dn, on=['Lon', 'Lat', 'Year'])
    d = d[(d.Year >= y0) & (d.Year < y0 + NYR)]
    g = d.groupby(['Lon', 'Lat'])
    first = d[d.Year == y0].set_index(['Lon', 'Lat'])
    last = d[d.Year == y0 + NYR - 1].set_index(['Lon', 'Lat'])
    out = pd.DataFrame({
        'bns0': first.BNS, 'bns19': last.BNS,
        'dens0': first.densBNS, 'dens19': last.densBNS,
        'agdd5': g.agdd5.mean(), 'grassPAR': g.grassPAR.mean(), 'mTmin20': g.mTmin20.mean(),
        'grass_mean': g.GRASSFPC.mean()}).dropna()
    out = out[out.bns0 >= 0.05]
    out['ret'] = out.bns19 / out.bns0
    out['ret_dens'] = out.dens19 / out.dens0
    return out.reset_index()


def binned(x, y, edges):
    idx = np.digitize(x, edges) - 1
    m = np.array([np.median(y[idx == i]) if (idx == i).sum() >= 8 else np.nan
                  for i in range(len(edges) - 1)])
    n = np.array([(idx == i).sum() for i in range(len(edges) - 1)])
    return 0.5 * (edges[1:] + edges[:-1]), m, n


def main():
    png = sys.argv[1]
    runs = {}
    for a in sys.argv[2:]:
        tag, rest = a.split('=', 1)
        root, y0 = rest.rsplit(':', 1)
        runs[tag] = cells(root, int(y0))
    edges = np.arange(300, 1500, 100.0)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    print(f'{"run":8s} {"cells":>6s} {"BNS ret":>8s} {"dens ret":>9s} {"agdd5":>7s} {"grassPAR":>10s}')
    for tag, c in runs.items():
        w = c.bns0
        print(f'{tag:8s} {len(c):6d} {np.average(c.ret, weights=w):8.3f} '
              f'{np.average(c.ret_dens, weights=w):9.3f} {c.agdd5.median():7.0f} {c.grassPAR.median():10.3g}')
        xm, ym, n = binned(c.agdd5.values, c.ret.values, edges)
        ax[0].scatter(c.agdd5, c.ret, s=6, alpha=0.25, color=COL.get(tag))
        ax[0].plot(xm, ym, '-o', lw=2, color=COL.get(tag), label=f'{tag} (n={len(c)})')
        ax[1].hist(c.agdd5, bins=edges, histtype='step', lw=2, color=COL.get(tag), label=tag)
        xm, ym, n = binned(c.grassPAR.values / 1e6, c.ret.values, np.arange(3.0, 8.0, 0.4))
        ax[2].plot(xm, ym, '-o', lw=2, color=COL.get(tag), label=tag)
        print('   by GDD5 bin (median retention, cells): ' + '  '.join(
            f'{int(e)}:{m:.2f}({k})' for e, m, k in zip(edges[:-1], binned(c.agdd5.values, c.ret.values, edges)[1],
                                                       binned(c.agdd5.values, c.ret.values, edges)[2]) if k >= 8))
    ax[0].set(xlabel='GDD5, 20-yr mean per cell', ylabel='BNS FPC year 19 / year 0',
              title='Larch retention against the cell\'s own warmth', ylim=(0, 1.6))
    ax[0].axhline(1, c='k', lw=0.6)
    ax[1].set(xlabel='GDD5', ylabel='cells', title='Where each run\'s larch cells sit')
    ax[2].set(xlabel='grassPAR from est_limits.out, 10$^6$ (LPJ-GUESS units), 20-yr mean', ylabel='median retention',
              title='Retention against light reaching the grass')
    for a in ax:
        a.legend(fontsize=8)
    fig.suptitle('Siberian box 55-75N 60-180E, cells with BNS FPC >= 0.05 at the start')
    fig.tight_layout()
    fig.savefig(png, dpi=130)


if __name__ == '__main__':
    main()
