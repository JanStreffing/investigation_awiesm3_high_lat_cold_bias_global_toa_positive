"""First decade of PICAL_crunveg_kpp_gm1000 against the TKE-only line and the IDEMIX arm, same years.

All three share the state of PICAL_crunveg_tke_albsn082 at the end of 2169 (the two arms start with a
cold atmosphere in 2170, which is discarded). kpp_gm1000 changes two things against the TKE-only line:
mix_scheme 'cvmix_TKE' -> 'KPP' and K_GM_max 2500 south / 1000 north -> 1000 everywhere.
INTERIM: the arm is still running; its years are read from the run's work directory through the view
reval/views/kpp_gm1000_interim.

Means over 2171-2179 with diff = arm minus TKE-only and 2 s.e. from the interannual spread; for the
temperature bands also the linear trend over the nine years. Labrador: March mean of MLD2 and ice
fraction in the interior box 56-62N 60-50W (labrador_any_mesh.py).

Usage: kpp_gm1000_first_decade.py [y0 y1]   (reval environment; writes data/clim/kpp_gm1000_first_decade.txt)
"""
import os, sys
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import pandas as pd
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
RUNS = {
    'TKE': dict(tag='tke_albsn082', root=f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082', lab='labrador_tke_albsn082_2160-2209.csv'),
    'IDEMIX': dict(tag='idemix_albsn082', root=f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_idemix_albsn082', lab='labrador_idemix_albsn082_2170-2209.csv'),
    'KPP+GM1000': dict(tag='kpp_gm1000', root=f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_kpp_gm1000', lab='labrador_kpp_gm1000_2170-2189.csv'),
}
Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (2171, 2179)
ACC, RHO_LF = 3600.0, 3.3355e8
BANDS = [(0, 100), (100, 300), (300, 700), (700, 1500), (1500, 2500), (2500, 6000)]
ZI = np.abs(xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')['nz'].values)
DZ, ZM = np.diff(ZI), 0.5 * (ZI[1:] + ZI[:-1])


def gmean(root, var, yr):
    with xr.open_dataset(f'{root}/outdata/oifs/atm_remapped_1m_{var}_{yr}-{yr}.nc') as ds:
        x = ds[var].mean([d for d in ds[var].dims if d not in ('lat', 'lon')]).mean('lon')
        w = np.cos(np.deg2rad(ds['lat']))
        return float((x * w).sum() / w.sum())


def table(r):
    s = pd.read_csv(f"{REPO}/data/clim/series_{r['tag']}.csv").set_index('year')
    lab = pd.read_csv(f"{REPO}/data/clim/{r['lab']}", header=None, names=['exp', 'year', 'month', 'region', 'var', 'val'])
    lab = lab[(lab.region == 'interior') & (lab.month == 3)].pivot(index='year', columns='var', values='val')
    rows = []
    for y in range(Y0, Y1 + 1):
        f = {v: gmean(r['root'], v, y) / ACC for v in ('tsr', 'ttr', 'ssr', 'str', 'sshf', 'slhf', 'sf')}
        prof = np.load(f"{P}/reval/spinup_cache/{r['tag']}/hovm_temp_{y}.npy")
        am = np.load(f"{P}/reval/spinup_cache/{r['tag']}/amoc_{y}.npy")
        row = dict(year=y, toa=f['tsr'] + f['ttr'], sfc=f['ssr'] + f['str'] + f['sshf'] + f['slhf'] - RHO_LF * f['sf'],
                   asr=f['tsr'], olr=-f['ttr'], amoc26=am[0], amoc_max=am[1],
                   lab_mld=abs(lab.loc[y, 'MLD2']), lab_ice=lab.loc[y, 'a_ice'])
        n = len(prof)
        for a, b in BANDS:
            m = (ZM[:n] >= a) & (ZM[:n] < b)
            row[f'T{a}-{b}'] = float(np.sum(prof[m] * DZ[:n][m]) / np.sum(DZ[:n][m]))
        rows.append(row)
    return pd.DataFrame(rows).set_index('year').join(s.drop(columns=['exp', 'tsr_raw', 'ttr_raw'], errors='ignore'))


ROWS = [('toa', 'net TOA', 'W/m2', 1), ('sfc', 'net surface', 'W/m2', 1), ('asr', 'absorbed solar', 'W/m2', 1),
        ('olr', 'outgoing longwave', 'W/m2', 1),
        ('t2m_glob', 'T2m global', 'C', 1), ('t2m_nh', 'T2m NH', 'C', 1), ('t2m_sh', 'T2m SH', 'C', 1),
        ('t2m_6090n', 'T2m 60-90N', 'C', 1), ('t2m_6090s', 'T2m 60-90S', 'C', 1),
        ('sivoln_m04', 'NH ice volume April', '1e3 km3', 1e-3), ('sivoln_m09', 'NH ice volume September', '1e3 km3', 1e-3),
        ('sivols_m09', 'SH ice volume September', '1e3 km3', 1e-3), ('sivols_m02', 'SH ice volume February', '1e3 km3', 1e-3),
        ('siextentn_m03', 'NH ice extent March', '1e6 km2', 1), ('siextents_m09', 'SH ice extent September', '1e6 km2', 1),
        ('siextents_m02', 'SH ice extent February', '1e6 km2', 1),
        ('amoc26', 'AMOC 26.5N', 'Sv', 1), ('amoc_max', 'AMOC max 20-60N', 'Sv', 1),
        ('lab_mld', 'Labrador March MLD', 'm', 1), ('lab_ice', 'Labrador March ice frac', '-', 1)] + \
       [(f'T{a}-{b}', f'T {a}-{b} m', 'C', 1) for a, b in BANDS]

T = {k: table(r) for k, r in RUNS.items()}
out = [f'PICAL_crunveg_kpp_gm1000 (KPP, uniform K_GM_max 1000) against the TKE-only line and the IDEMIX arm, {Y0}-{Y1}.',
       'Means; d = arm minus TKE-only, +- 2 s.e. of the difference; * marks |d| > 2 s.e.',
       f"{'quantity':26s} {'unit':8s} {'TKE':>9s} {'IDEMIX':>9s} {'KPP+GM1000':>11s} {'d IDEMIX':>10s} {'d KPP+GM1000':>18s}"]
for key, name, unit, sc in ROWS:
    a = T['TKE'][key].values * sc
    line = f"{name:26s} {unit:8s} {np.nanmean(a):9.3f}"
    ds = []
    for k in ('IDEMIX', 'KPP+GM1000'):
        b = T[k][key].values * sc
        d = np.nanmean(b) - np.nanmean(a)
        se = 2 * np.sqrt(np.nanvar(a, ddof=1) / len(a) + np.nanvar(b, ddof=1) / len(b))
        line += f" {np.nanmean(b):{9 if k == 'IDEMIX' else 11}.3f}"
        ds.append(f"{d:+8.3f}{'*' if abs(d) > se else ' '}" + (f" +-{se:6.3f}" if k != 'IDEMIX' else ''))
    out.append(line + '  ' + '  '.join(ds))
out += ['', f'Temperature trend {Y0}-{Y1}, mK/decade +- 2 s.e.', f"{'band':14s} {'TKE':>16s} {'IDEMIX':>16s} {'KPP+GM1000':>16s}"]
yy = np.arange(Y0, Y1 + 1, dtype=float)
for a, b in BANDS:
    line = f"{f'{a}-{b} m':14s}"
    for k in T:
        x = T[k][f'T{a}-{b}'].values
        A = np.vstack([yy - yy.mean(), np.ones_like(yy)]).T
        c, res = np.linalg.lstsq(A, x, rcond=None)[:2]
        se = np.sqrt(res[0] / (len(yy) - 2) / np.sum((yy - yy.mean()) ** 2))
        line += f" {1e4 * c[0]:+8.1f} +-{2e4 * se:5.1f}"
    out.append(line)
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/kpp_gm1000_{Y0}-{Y1}.txt', 'w').write(txt + '\n')
