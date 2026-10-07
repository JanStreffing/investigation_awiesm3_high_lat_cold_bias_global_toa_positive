"""Sea-ice volume and ocean drift of the four arms of 2026-10-06/07 against PICAL_crunveg_tke_albsn082.

Arms (each one change from the control at the 2169 restart): icb (Antarctic calving as icebergs),
wamz0 (wave-model roughness on the ECWAM default), cdoi10 (ice-ocean drag 0.0055 -> 0.010),
cpres15 (ice strength c_pressure 20 -> 15). 2170 is the arms' atmosphere spin-up year and is left out.

Decadal means 2171-2179 and 2180-2189 of sea-ice volume and extent, T2m, the energy budget and the
Atlantic overturning; arm minus control, with 2 standard errors of the difference from the interannual
spread of both runs (no autocorrelation correction; * marks a difference beyond it). Then the linear
trend 2171-2189 of the global-mean ocean temperature by depth band, mK per decade +- 2 s.e., and the
band means of 2180-2189 minus the control.

Usage: overnight_arms_volumes_drift.py [runs periods outname]
  runs     comma list of control plus arms, default control,icb,wamz0,cdoi10,cpres15; also ctl2 (second
           realisation of the control), z0def, z0floor, z0form (sea-ice roughness runs; their control is z0def:
           the first name in the list is the reference)
  periods  comma list like 2171-2179,2180-2189; the trends use the full span of the periods
  outname  data/clim/<outname>.txt, default overnight_arms_volumes_drift
(reval environment; needs series_<tag>.csv and the spinup cache of each run)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import pandas as pd
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
import sys
ALL = {'control': ('tke_albsn082', 'PICAL_crunveg_tke_albsn082'), **{k: (k, f'PICAL_crunveg_{k}') for k in
       ('icb', 'wamz0', 'cdoi10', 'cpres15', 'ctl2', 'z0def', 'z0floor', 'z0form')}}
SEL = sys.argv[1].split(',') if len(sys.argv) > 1 else ['control', 'icb', 'wamz0', 'cdoi10', 'cpres15']
RUNS = [(k, *ALL[k]) for k in SEL]
REF = SEL[0]
PER = [tuple(int(v) for v in p.split('-')) for p in sys.argv[2].split(',')] if len(sys.argv) > 2 else [(2171, 2179), (2180, 2189)]
OUT = sys.argv[3] if len(sys.argv) > 3 else 'overnight_arms_volumes_drift'
Y0, Y1 = PER[0][0], PER[-1][1]
ACC, RHO_LF = 3600.0, 3.3355e8
BANDS = [(0, 100), (100, 300), (300, 700), (700, 1500), (1500, 2500), (2500, 6000)]
ZI = np.abs(xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')['nz'].values)
DZ, ZM = np.diff(ZI), 0.5 * (ZI[1:] + ZI[:-1])


def gmean(root, var, yr):
    with xr.open_dataset(f'{root}/outdata/oifs/atm_remapped_1m_{var}_{yr}-{yr}.nc') as ds:
        x = ds[var].mean([d for d in ds[var].dims if d not in ('lat', 'lon')]).mean('lon')
        w = np.cos(np.deg2rad(ds['lat']))
        return float((x * w).sum() / w.sum())


def table(tag, exp):
    root = f'{P}/runtime/awiesm3-v3.4/{exp}'
    s = pd.read_csv(f'{REPO}/data/clim/series_{tag}.csv').set_index('year')
    rows = []
    for y in range(Y0, Y1 + 1):
        f = {v: gmean(root, v, y) / ACC for v in ('tsr', 'ttr', 'ssr', 'str', 'sshf', 'slhf', 'sf')}
        prof = np.load(f'{P}/reval/spinup_cache/{tag}/hovm_temp_{y}.npy')
        am = np.load(f'{P}/reval/spinup_cache/{tag}/amoc_{y}.npy')
        row = dict(year=y, toa=f['tsr'] + f['ttr'], sfc=f['ssr'] + f['str'] + f['sshf'] + f['slhf'] - RHO_LF * f['sf'],
                   amoc26=am[0], amoc_max=am[1])
        n = len(prof)
        for a, b in BANDS:
            m = (ZM[:n] >= a) & (ZM[:n] < b)
            row[f'T{a}-{b}'] = float(np.sum(prof[m] * DZ[:n][m]) / np.sum(DZ[:n][m]))
        rows.append(row)
    return pd.DataFrame(rows).set_index('year').join(s.drop(columns=['exp', 'tsr_raw', 'ttr_raw'], errors='ignore'))


ROWS = [('sivoln_m04', 'NH ice volume April', '1e3 km3', 1e-3), ('sivoln_m09', 'NH ice volume September', '1e3 km3', 1e-3),
        ('sivols_m09', 'SH ice volume September', '1e3 km3', 1e-3), ('sivols_m02', 'SH ice volume February', '1e3 km3', 1e-3),
        ('siextentn_m03', 'NH ice extent March', '1e6 km2', 1), ('siextentn_m09', 'NH ice extent September', '1e6 km2', 1),
        ('siextents_m09', 'SH ice extent September', '1e6 km2', 1), ('siextents_m02', 'SH ice extent February', '1e6 km2', 1),
        ('t2m_glob', 'T2m global', 'C', 1), ('t2m_6090n', 'T2m 60-90N', 'C', 1), ('t2m_6090s', 'T2m 60-90S', 'C', 1),
        ('toa', 'net TOA', 'W/m2', 1), ('sfc', 'net surface', 'W/m2', 1),
        ('amoc26', 'AMOC 26.5N', 'Sv', 1), ('amoc_max', 'AMOC max 20-60N', 'Sv', 1)]
T = {k: table(tag, exp) for k, tag, exp in RUNS}
arms = [k for k, _, _ in RUNS if k != REF]
out = [f'Runs from the 2169 restart against {ALL[REF][1]} ({REF}). Period means; the other columns are run minus {REF};',
       '* marks a difference beyond 2 standard errors of the difference (interannual spread of both runs).', '']
for p0, p1 in PER:
    out.append(f'{p0}-{p1}')
    out.append(f"{'quantity':26s}{'unit':9s}{REF:>9s}" + ''.join(f'{a:>11s}' for a in arms))
    for key, name, unit, sc in ROWS:
        c = T[REF].loc[p0:p1, key].values * sc
        line = f'{name:26s}{unit:9s}{np.nanmean(c):9.3f}'
        for a in arms:
            x = T[a].loc[p0:p1, key].values * sc
            d = np.nanmean(x) - np.nanmean(c)
            se = 2 * np.sqrt(np.nanvar(c, ddof=1) / len(c) + np.nanvar(x, ddof=1) / len(x))
            line += f"{d:+10.3f}{'*' if abs(d) > se else ' '}"
        out.append(line)
    out.append('')
yy = np.arange(Y0, Y1 + 1, dtype=float)
out.append(f'Global-mean ocean temperature: trend {Y0}-{Y1}, mK/decade (+- 2 s.e.)')
out.append(f"{'band':14s}" + ''.join(f'{k:>16s}' for k, _, _ in RUNS))
for a, b in BANDS:
    line = f"{f'{a}-{b} m':14s}"
    for k, _, _ in RUNS:
        x = T[k][f'T{a}-{b}'].values
        A = np.vstack([yy - yy.mean(), np.ones_like(yy)]).T
        c, res = np.linalg.lstsq(A, x, rcond=None)[:2]
        se = np.sqrt(res[0] / (len(yy) - 2) / np.sum((yy - yy.mean()) ** 2))
        line += f'{1e4 * c[0]:+9.1f}+-{2e4 * se:5.1f}'
    out.append(line)
out += ['', f'Band mean {PER[-1][0]}-{PER[-1][1]}, run minus {REF}, mK', f"{'band':14s}" + ''.join(f'{a:>11s}' for a in arms)]
for a, b in BANDS:
    out.append(f"{f'{a}-{b} m':14s}" + ''.join(f"{1e3 * (T[k].loc[PER[-1][0]:PER[-1][1], f'T{a}-{b}'].mean() - T[REF].loc[PER[-1][0]:PER[-1][1], f'T{a}-{b}'].mean()):+11.1f}" for k in arms))
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/{OUT}.txt', 'w').write(txt + '\n')
