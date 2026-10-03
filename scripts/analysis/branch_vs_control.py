"""Branches against their control over the same years: net TOA, net surface energy, T2m.

For each run, per year (12 complete months of OpenIFS remapped output, from outdata or a
running leg's work directory): global net TOA, net energy into the surface
(ssr+str+slhf+sshf minus the latent heat of snowfall, sf*3.3355e8), and T2m globally, by
band, and 60-90S in JJA.  Differences against CONTROL are taken over the years every run
has.  The detection threshold of a difference of N-year means is 2*sd*sqrt(2/N), with sd
the control's year-to-year standard deviation over those years; it ignores autocorrelation,
so it is a floor.  Branches that cold-start the atmosphere are also compared without their
first year (SKIP1).

Usage:  CONTROL=PI200 RUNS=PI200_spp,PI200_mle,PI200_h0,PI200_all3 Y0=1560 Y1=1579 \
        python3 scripts/analysis/branch_vs_control.py
"""
import os, glob
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
CONTROL = os.environ.get('CONTROL', 'PI200')
RUNS = os.environ.get('RUNS', 'PI200_spp,PI200_mle,PI200_h0,PI200_all3').split(',')
Y0, Y1 = int(os.environ.get('Y0', 1560)), int(os.environ.get('Y1', 1579))
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'
BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45), ('30S-30N', -30, 30), ('45-60N', 45, 60), ('60-90N', 60, 90)]


def path(arm, var, y):
    work = [q for q in sorted(glob.glob(f'{R}/{arm}/run_*/work/atm_remapped_1m_{var}_{y}-{y}.nc'))
            if os.path.basename(os.path.dirname(os.path.dirname(q))).count('.') == 0]
    for p in [f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'] + work:
        if os.path.exists(p):
            return p
    return None


with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = lsm[0] if lsm.ndim == 3 else lsm
    lat = np.squeeze(d['lat'].values)
LAT = np.broadcast_to(lat[:, None], lsm.shape); WT = np.cos(np.deg2rad(LAT))


def g(a, k=None):
    k = np.ones(a.shape, bool) if k is None else k
    return float(np.average(a[k], weights=WT[k]))


def year(arm, y):
    f = {}
    for v in ('tsr', 'ttr', 'ssr', 'str', 'slhf', 'sshf', 'sf', '2t'):
        p = path(arm, v, y)
        if p is None:
            return None
        with xr.open_dataset(p, decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8')
        if a.shape[0] != 12 or not np.isfinite(a[-1]).any():
            return None
        f[v] = a
    ann = {v: f[v].mean(0) for v in f}
    for v in ('tsr', 'ttr', 'ssr', 'str', 'slhf', 'sshf'):
        ann[v] = ann[v] / 3600.0
    out = {'net TOA': g(ann['tsr'] + ann['ttr']),
           'sfc net': g(ann['ssr'] + ann['str'] + ann['slhf'] + ann['sshf'] - ann['sf'] * 1000 * 333550.0 / 3600.0),
           'T2m glob': g(ann['2t'])}
    for n, a, b in BANDS:
        out[f'T2m {n}'] = g(ann['2t'], (LAT >= a) & (LAT < b))
    jja = f['2t'][[5, 6, 7]].mean(0)
    out['T2m 90-60S JJA'] = g(jja, LAT < -60)
    out['T2m 90-60S JJA ocn'] = g(jja, (LAT < -60) & (lsm <= 0.5))
    return out


D = {}
for arm in [CONTROL] + RUNS:
    D[arm] = {}
    for y in range(Y0, Y1 + 1):
        r = year(arm, y)
        if r is not None:
            D[arm][y] = r
    print(f'{arm}: complete years {min(D[arm]) if D[arm] else "-"}-{max(D[arm]) if D[arm] else "-"} ({len(D[arm])})')

KEYS = list(next(iter(D[CONTROL].values())).keys())
common = sorted(set.intersection(*[set(D[a]) for a in D]))
print(f'\ncommon years {common[0]}-{common[-1]} (N={len(common)})')

print('\nper-year global net TOA / sfc net / T2m glob / T2m 90-60S JJA ocean')
print(f'  {"year":<6}' + ''.join(f'{a:>28}' for a in D))
for y in sorted(set().union(*[set(D[a]) for a in D])):
    row = f'  {y:<6}'
    for a in D:
        r = D[a].get(y)
        row += f'{"":>28}' if r is None else f'{r["net TOA"]:+7.2f}{r["sfc net"]:+7.2f}{r["T2m glob"]:7.2f}{r["T2m 90-60S JJA ocn"]:7.2f}'
    print(row)

for label, yrs in (('all common years', common), ('SKIP1: without the first year', common[1:])):
    N = len(yrs)
    ctl = {k: np.array([D[CONTROL][y][k] for y in yrs]) for k in KEYS}
    print(f'\n{label}: {yrs[0]}-{yrs[-1]}, N={N}.  Control mean, branch minus control, threshold 2*sd*sqrt(2/N)')
    print(f'  {"":<20}{CONTROL:>10}' + ''.join(f'{a.replace(CONTROL + "_", ""):>9}' for a in RUNS) + f'{"thresh":>9}')
    for k in KEYS:
        c = ctl[k].mean(); thr = 2 * ctl[k].std(ddof=1) * np.sqrt(2.0 / N)
        cells = ''
        for a in RUNS:
            dlt = np.mean([D[a][y][k] for y in yrs]) - c
            cells += f'{dlt:+8.2f}{"*" if abs(dlt) > thr else " "}'
        print(f'  {k:<20}{c:10.2f}{cells}{thr:9.2f}')
print('\n* beyond the threshold.  T2m in K, fluxes in W/m2 (positive down).')
