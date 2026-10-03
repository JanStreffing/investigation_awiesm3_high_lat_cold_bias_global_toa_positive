"""Where does GM 2500 take heat uptake away? Pairs that share a start state:
15F (GM 1000) / 16C (1500) / 16E (2500), all from the 1350 restart, and
PI200_gmR1500 / gmR2500, both from PI200 at 1579-12-31.
Heat uptake per region and depth band = OHC change over the common period / time,
in W per m2 of Earth surface (so regions add up to the global number). Surface flux the same,
uptake = -fh. Input: data/clim/ohc/*.csv from ohc_by_region_depth.py.
"""
import glob, numpy as np, pandas as pd
d = pd.concat(pd.read_csv(f, names=['exp', 'year', 'kind', 'region', 'band', 'v']) for f in glob.glob('data/clim/ohc/*.csv'))
AE, YR = 5.1e14, 365.25 * 86400
REG = ['SO', 'SPNA', 'NORDARC', 'REST']; BAND = ['0-700', '700-2000', '2000+']
def uptake(e, a, b):
    o = d[(d.exp == e) & (d.kind == 'ohc')].pivot_table(index='year', columns=['region', 'band'], values='v')
    # annual means are centred mid-year: change from year a to year b spans b-a years
    r = (o.loc[b] - o.loc[a]) / ((b - a) * YR) / AE
    f = -d[(d.exp == e) & (d.kind == 'fh') & (d.year > a) & (d.year <= b)].groupby('region').v.mean() / AE
    return r, f
def table(title, arms, a, b):
    print(f'\n=== {title}: OHC change {a} -> {b}, W/m2 of Earth surface ===')
    U = {e: uptake(e, a, b) for e in arms}
    hdr = 'region  band      ' + ''.join(f'{e:>15s}' for e in arms)
    for i in range(1, len(arms)): hdr += f'{arms[i]+"-"+arms[0]:>26s}'
    print(hdr)
    for r in REG:
        for bnd in BAND:
            vals = [U[e][0][(r, bnd)] for e in arms]
            print(f'{r:7s} {bnd:9s} ' + ''.join(f'{v:15.3f}' for v in vals) + ''.join(f'{vals[i]-vals[0]:26.3f}' for i in range(1, len(arms))))
    for lab, sel in (('TOTAL', None), ('SO all', 'SO'), ('SPNA all', 'SPNA'), ('NORDARC all', 'NORDARC'), ('REST all', 'REST')):
        vals = [U[e][0].sum() if sel is None else U[e][0][sel].sum() for e in arms]
        print(f'{lab:17s} ' + ''.join(f'{v:15.3f}' for v in vals) + ''.join(f'{vals[i]-vals[0]:26.3f}' for i in range(1, len(arms))))
    print('surface uptake (-fh), period mean:')
    for r in REG + ['all']:
        vals = [U[e][1].sum() if r == 'all' else U[e][1][r] for e in arms]
        print(f'  {r:15s} ' + ''.join(f'{v:15.3f}' for v in vals) + ''.join(f'{vals[i]-vals[0]:26.3f}' for i in range(1, len(arms))))
table('GM dose, shared 1350 start, first 30 years', ['15F', '16C', '16E'], 1350, 1379)
table('GM 1500 vs 2500, shared 1350 start, 40 years', ['16C', '16E'], 1350, 1389)
table('Rossby-scaled GM 1500 vs 2500 from PI200 1579', ['PI200_gmR1500', 'PI200_gmR2500'], 1580, 1599)
