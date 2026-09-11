"""Arm-minus-control differences for clean pairs, from the two annual stores.

For each pair (control, arm) over their common complete years: heating rate by depth range
and in total (W/m2 of Earth, from annual heat content), net TOA, and the AMOC indices, as
5-yr chunk means and as the mean over all common years with the standard error from the
year-to-year scatter of the annual difference.

Usage:  PAIRS=15F:16B,16C:16D python3 scripts/analysis/pair_diff_store.py
"""
import os
import numpy as np, xarray as xr
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AE, SY = 5.101e14, 365.25 * 86400
C = xr.open_dataset(os.path.join(REPO, 'data', 'coupled_annual_diag.nc'))
A = xr.open_dataset(os.path.join(REPO, 'data', 'amoc_annual_diag.nc'))
PAIRS = [p.split(':') for p in os.environ.get('PAIRS', '15F:16B,16C:16D').split(',')]
DEP = [('0-100', 'ohc_0_100'), ('100-700', 'ohc_100_700'), ('700-2000', 'ohc_700_2000'), ('>2000', 'ohc_gt2000')]

def annual_rates(run):
    s = C.sel(run=run); y = C['year'].values.astype(int)
    H = {n: s[v].values for n, v in DEP}; H['total'] = sum(H[n] for n, _ in DEP)
    out = {n: dict(zip(y[1:], np.diff(h) / (SY * AE))) for n, h in H.items()}   # rate over year y-1 -> y
    out['net TOA'] = dict(zip(y, s['toa'].values))
    if run in A['run'].values.astype(str):
        a = A.sel(run=run); ya = A['year'].values.astype(int)
        for v in ('amoc26', 'amoc26_res', 'amoc4060', 'amoc4060_res'):
            out[v] = dict(zip(ya, a[v].values))
    return out

for ctl, arm in PAIRS:
    rc, ra = annual_rates(ctl), annual_rates(arm)
    print(f'\n=== {arm} minus {ctl}')
    for n in ['0-100', '100-700', '700-2000', '>2000', 'total', 'net TOA', 'amoc26', 'amoc26_res', 'amoc4060', 'amoc4060_res']:
        if n not in rc or n not in ra:
            continue
        ys = sorted(y for y in set(rc[n]) & set(ra[n]) if np.isfinite(rc[n][y]) and np.isfinite(ra[n][y]))
        if len(ys) < 2:
            continue
        d = np.array([ra[n][y] - rc[n][y] for y in ys]); y0 = ys[0]
        chunks = []
        for k in range(0, 50, 5):
            sel = [i for i, y in enumerate(ys) if y0 + k <= y < y0 + k + 5]
            if sel:
                chunks.append(f'{d[sel].mean():+6.2f}' + ('' if len(sel) == 5 else f'({len(sel)})'))
        print(f'  {n:<13}{ys[0]}-{ys[-1]}  mean {d.mean():+.3f} se {d.std(ddof=1)/np.sqrt(len(d)):.3f}   5-yr: ' + ' '.join(chunks))
