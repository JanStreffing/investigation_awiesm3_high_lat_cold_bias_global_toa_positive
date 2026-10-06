"""Validation of the v3.5.0 candidate build against the run it repeats, over the same years.

PICAL_v350val repeats 2190-2199 of PICAL_crunveg_tke_albsn082 from the same restart with every
component built from the tagged commits (notes/awiesm3_v3.5.0_component_tags.md).  The two runs
share a start and diverge chaotically, so a year of one is not the same weather as that year of
the other: the comparison is between the two 9-year means (2191-2199, 2190 is the cold-started
atmosphere of the candidate run), with 2 standard errors of the difference from the interannual
spread of each run (no autocorrelation correction).

Quantities as in tke_vs_idemix_branch_score.py: TOA and surface budgets from the remapped monthly
OpenIFS files, T2m and sea ice from series_<run>.csv, temperature by depth band and AMOC from the
spin-up cache.

Usage: v350_validation_pairs.py     (reval environment; writes data/clim/v350_validation_pairs.txt)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import pandas as pd
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
RUNS = {
    'tke_albsn082': dict(root=f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082', cache=f'{P}/reval/spinup_cache/tke_albsn082'),
    'v350val': dict(root=f'{P}/runtime/awiesm3-v3.4/PICAL_v350val', cache=f'{P}/reval/spinup_cache/v350val'),
}
Y0, Y1 = 2191, 2199
ACC, RHO_LF = 3600.0, 3.3355e8
BANDS = [(0, 100), (100, 300), (300, 700), (700, 1500), (1500, 2500), (2500, 6000)]
ZI = np.abs(xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')['nz'].values)
DZ, ZM = np.diff(ZI), 0.5 * (ZI[1:] + ZI[:-1])


def gmean(root, var, yr):
    with xr.open_dataset(f'{root}/outdata/oifs/atm_remapped_1m_{var}_{yr}-{yr}.nc') as ds:
        x = ds[var].mean([d for d in ds[var].dims if d not in ('lat', 'lon')])
        w = np.cos(np.deg2rad(ds['lat']))
        return float((x * w).sum(('lat', 'lon')) / (w.sum('lat') * ds.sizes['lon']))


def table(tag):
    r = RUNS[tag]
    s = pd.read_csv(f'{REPO}/data/clim/series_{tag}.csv').set_index('year')
    rows = []
    for y in range(Y0, Y1 + 1):
        f = {v: gmean(r['root'], v, y) / ACC for v in ('tsr', 'ttr', 'ssr', 'str', 'sshf', 'slhf', 'sf')}
        prof = np.load(f"{r['cache']}/hovm_temp_{y}.npy")
        am = np.load(f"{r['cache']}/amoc_{y}.npy")
        row = dict(year=y, toa=f['tsr'] + f['ttr'], sfc=f['ssr'] + f['str'] + f['sshf'] + f['slhf'] - RHO_LF * f['sf'],
                   asr=f['tsr'], olr=-f['ttr'], amoc26=am[0], amoc_max=am[1])
        n = len(prof)
        for a, b in BANDS:
            m = (ZM[:n] >= a) & (ZM[:n] < b)
            row[f'T{a}-{b}'] = float(np.sum(prof[m] * DZ[:n][m]) / np.sum(DZ[:n][m]))
        rows.append(row)
    d = pd.DataFrame(rows).set_index('year')
    return d.join(s.drop(columns=['exp', 'tsr_raw', 'ttr_raw'], errors='ignore'))


ROWS = [('toa', 'net TOA', 'W/m2', 1), ('sfc', 'net surface', 'W/m2', 1), ('asr', 'absorbed solar', 'W/m2', 1),
        ('olr', 'outgoing longwave', 'W/m2', 1),
        ('t2m_glob', 'T2m global', 'C', 1), ('t2m_nh', 'T2m NH', 'C', 1), ('t2m_sh', 'T2m SH', 'C', 1),
        ('t2m_6090n', 'T2m 60-90N', 'C', 1), ('t2m_6090s', 'T2m 60-90S', 'C', 1),
        ('t2m_land', 'T2m land', 'C', 1), ('t2m_ocean', 'T2m ocean', 'C', 1),
        ('sivoln_m04', 'NH ice volume April', '1e3 km3', 1e-3), ('sivoln_m09', 'NH ice volume September', '1e3 km3', 1e-3),
        ('sivols_m09', 'SH ice volume September', '1e3 km3', 1e-3), ('sivols_m02', 'SH ice volume February', '1e3 km3', 1e-3),
        ('siextentn_m03', 'NH ice extent March', '1e6 km2', 1), ('siextentn_m09', 'NH ice extent September', '1e6 km2', 1),
        ('siextents_m09', 'SH ice extent September', '1e6 km2', 1), ('siextents_m02', 'SH ice extent February', '1e6 km2', 1),
        ('amoc26', 'AMOC 26.5N', 'Sv', 1), ('amoc_max', 'AMOC max 20-60N', 'Sv', 1)] + \
       [(f'T{a}-{b}', f'T {a}-{b} m', 'C', 1) for a, b in BANDS]

A, B = table('tke_albsn082'), table('v350val')
out = [f'PICAL_v350val (v3.5.0 candidate build) against PICAL_crunveg_tke_albsn082 (tuning build), same restart, {Y0}-{Y1}.',
       'Means of the 9 years; diff = candidate minus tuning build; +- is 2 standard errors of the difference; * marks |diff| > 2 s.e.',
       f"{'quantity':30s} {'unit':8s} {'tuning':>10s} {'candidate':>10s} {'diff':>9s} {'+- 2 s.e.':>10s}"]
for key, name, unit, sc in ROWS:
    if key not in A or key not in B:
        continue
    a, b = A[key].values * sc, B[key].values * sc
    d = np.nanmean(b) - np.nanmean(a)
    se = 2 * np.sqrt(np.nanvar(a, ddof=1) / len(a) + np.nanvar(b, ddof=1) / len(b))
    out.append(f"{name:30s} {unit:8s} {np.nanmean(a):10.3f} {np.nanmean(b):10.3f} {d:+9.3f} {se:10.3f} {'*' if abs(d) > se else ''}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/v350_validation_pairs.txt', 'w').write(txt + '\n')
pd.concat([A.assign(run='tke_albsn082'), B.assign(run='v350val')]).to_csv(f'{REPO}/data/clim/v350_validation_pairs_annual.csv', float_format='%.6g')
