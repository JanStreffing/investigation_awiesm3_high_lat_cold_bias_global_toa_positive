"""Compare heat pathways between GM 1500 and GM 2500 for the three pairs, from the output of
heat_pathways.py (data/clim/heatpath/*.csv). Every difference is GM 1500 minus GM 2500:
  16C - 16E                       (1350 restart, 1351-1389)
  PI200_gmR1500 - PI200_gmR2500   (PI200 1580, Rossby scaling, 1581-1599)
  gm1500 - gm2500                 (crunveg 2119, 2121-2129; gm2500 = PICAL_crunveg_ob1200)
Per region: surface entry, lateral import (tend - sfc), storage by depth band (dOHC/dt);
Atlantic northward heat transport by latitude and depth band (resolved + bolus).
Units W m-2 of Earth surface (1 PW = 1.96 W m-2). '*' = beyond 2 sigma/sqrt(N) of the paired
annual difference series.
"""
import glob, numpy as np, pandas as pd
AE, YR = 5.1e14, 365.25 * 86400
d = pd.concat(pd.read_csv(f, names=['exp', 'year', 'kind', 'region', 'band', 'v']) for f in glob.glob('data/clim/heatpath/*.csv'))
d = d.drop_duplicates(['exp', 'year', 'kind', 'region', 'band'], keep='last')
PAIRS = [('16C', '16E', 1351, 1389), ('PI200_gmR1500', 'PI200_gmR2500', 1581, 1599), ('gm1500', 'gm2500', 2121, 2129)]
REG = ['SO', 'SATL', 'NATL', 'SPNA', 'NORD', 'IP', 'OTHER']
BANDS = ['0-700', '700-2000', '2000+']

def series(exp):
    x = d[d.exp == exp]
    o = x[x.kind == 'ohc'].pivot_table(index='year', columns=['region', 'band'], values='v').sort_index()
    rate = o.diff() / YR / AE                                             # storage, W m-2 Earth
    c = x[x.band == 'col'].pivot_table(index='year', columns=['kind', 'region'], values='v') / AE
    sfc, lat = c['sfc'], c['tend'] - c['sfc']
    m = x[x.kind.str.startswith('mht')].pivot_table(index='year', columns=['region', 'band'], values='v', aggfunc='sum') / AE
    return rate, sfc, lat, m

def diff(a, b, y0, y1):
    a, b = a.loc[y0:y1], b.loc[y0:y1]
    dd = (a - b).dropna(how='all')
    return dd.mean(), 2 * dd.std() / np.sqrt(dd.count())

def fmt(v, t): return f'{v:+8.3f}{"*" if abs(v) > t else " "}'

out = []
for ea, eb, y0, y1 in PAIRS:
    if ea not in set(d.exp) or eb not in set(d.exp):
        out.append(f'\n### {ea} - {eb}: missing data'); continue
    ra, sa, la, ma = series(ea); rb, sb, lb, mb = series(eb)
    out.append(f'\n### {ea} - {eb}  ({y0}-{y1}), W m-2 of Earth surface, GM 1500 minus GM 2500')
    out.append(f'{"region":7s} {"surface in":>10s} {"lateral in":>10s} {"store 0-700":>11s} {"700-2000":>9s} {"2000+":>9s} {"store all":>10s}')
    ms, ts = diff(sa, sb, y0, y1); ml, tl = diff(la, lb, y0, y1); mr, tr = diff(ra, rb, y0, y1)
    tot = np.zeros(6)
    for r in REG:
        st = [mr.get((r, b), np.nan) for b in BANDS] + [mr.get((r, 'all'), np.nan)]
        stt = [tr.get((r, b), np.nan) for b in BANDS] + [tr.get((r, 'all'), np.nan)]
        out.append(f'{r:7s} {fmt(ms[r], ts[r]):>10s} {fmt(ml[r], tl[r]):>10s} ' + ' '.join(f'{fmt(v, t):>9s}' for v, t in zip(st, stt)))
        tot += np.array([ms[r], ml[r]] + st)
    out.append(f'{"sum":7s} {tot[0]:+8.3f}   {tot[1]:+8.3f}   ' + ' '.join(f'{v:+8.3f} ' for v in tot[2:]))
    mm, tm = diff(ma, mb, y0, y1)
    lats = sorted({k[0] for k in mm.index}, key=lambda s: int(s.split('@')[1]))
    out.append('Atlantic northward heat transport (resolved+bolus), difference:')
    out.append(f'{"":9s}' + ''.join(f'{b:>11s}' for b in BANDS + ['all']))
    for la_ in lats:
        out.append(f'{la_:9s}' + ''.join(f'{fmt(mm[(la_, b)], tm[(la_, b)]):>11s}' for b in BANDS + ['all']))
    ref = mb.loc[y0:y1].mean()
    out.append(f'  (GM 2500 level, all depths: ' + ', '.join(f'{la_.split("@")[1]}: {ref[(la_, "all")]*AE/1e15:.2f} PW' for la_ in lats) + ')')
txt = '\n'.join(out); print(txt)
open('data/clim/heat_pathways_compare.txt', 'w').write(txt + '\n')
