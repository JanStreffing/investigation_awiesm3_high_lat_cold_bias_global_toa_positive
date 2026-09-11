"""Does the stable-BL lever RSBLB = 2 change ocean heat uptake by depth?

Two clean pairs in data/coupled_annual_diag.nc differ only in RSBLB (5 -> 2):
11N -> 11Q (1850 forcing) and 11P -> 11R (1990 forcing), 40 yr each on core3_beta.
Prints, per pair, the heating rate per depth range for each decade and over years 1-40
(W/m2 of Earth, from annual-mean heat content), decadal net TOA, and RSBLB-minus-control.
Year-to-year scatter of the annual rate difference gives the standard error of each mean.

Usage:  python3 scripts/analysis/rsblb_pair_ohc.py
"""
import os
import numpy as np, xarray as xr
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AE, SY = 5.101e14, 365.25 * 86400
ds = xr.open_dataset(os.path.join(REPO, 'data', 'coupled_annual_diag.nc'))
yrs = ds['year'].values.astype(int)
DEP = [('0-100', 'ohc_0_100'), ('100-700', 'ohc_100_700'), ('700-2000', 'ohc_700_2000'), ('>2000', 'ohc_gt2000')]

def series(run):
    s = ds.sel(run=run)
    H = {n: s[v].values for n, v in DEP}
    H['all'] = sum(H[n] for n, _ in DEP)
    return H, s['toa'].values

def rate(h, y0, y1):          # mean heating rate between annual means of y0 and y1
    i0, i1 = np.where(yrs == y0)[0][0], np.where(yrs == y1)[0][0]
    return (h[i1] - h[i0]) / ((y1 - y0) * SY * AE)

for ctl, arm, forcing in (('11N', '11Q', '1850'), ('11P', '11R', '1990')):
    Hc, tc = series(ctl); Ha, ta = series(arm)
    print(f'\n=== {ctl} -> {arm}  (RSBLB 5 -> 2, {forcing} forcing)   [W/m2 of Earth]')
    print(f'  {"":<10}' + ''.join(f'{d:>22}' for d in ('yrs 1-10', '11-20', '21-30', '31-40', '1-40')))
    spans = [(1350, 1360), (1360, 1370), (1370, 1380), (1380, 1389), (1350, 1389)]
    for n in [d[0] for d in DEP] + ['all']:
        cells = []
        for a, b in spans:
            rc, ra = rate(Hc[n], a, b), rate(Ha[n], a, b)
            cells.append(f'{rc:6.2f}{ra:6.2f} {ra - rc:+6.2f}')
        print(f'  {n:<10}' + ''.join(f'{c:>22}' for c in cells))
    cells = []
    for a, b in spans:
        k = (yrs >= a) & (yrs < b if b != 1389 else yrs <= b)
        cells.append(f'{np.nanmean(tc[k]):6.2f}{np.nanmean(ta[k]):6.2f} {np.nanmean(ta[k]) - np.nanmean(tc[k]):+6.2f}')
    print(f'  {"net TOA":<10}' + ''.join(f'{c:>22}' for c in cells))
    # standard error of the 40-yr mean difference from year-to-year scatter
    for n in ('700-2000', '>2000', 'all'):
        dc, da = np.diff(Hc[n]) / (SY * AE), np.diff(Ha[n]) / (SY * AE)
        k = np.isfinite(dc) & np.isfinite(da); dd = (da - dc)[k]
        print(f'  {n:<10} annual-rate difference over {k.sum()} yr: mean {dd.mean():+.3f}, se {dd.std(ddof=1)/np.sqrt(k.sum()):.3f}')
    dt = (ta - tc)[np.isfinite(ta - tc)]
    print(f'  net TOA    annual difference over {dt.size} yr: mean {dt.mean():+.3f}, se {dt.std(ddof=1)/np.sqrt(dt.size):.3f}')
print('\ncolumns per cell: control, RSBLB arm, difference')
