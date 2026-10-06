"""IDEMIX control arm against the TKE-only line, over the same years.

PICAL_crunveg_idemix_albsn082 branches from PICAL_crunveg_tke_albsn082 at the end of 2169 and
differs in one setting: mix_scheme 'cvmix_TKE+cvmix_IDEMIX' instead of 'cvmix_TKE'.  Both carry
albsn 0.82 and run to 2209.  The question is whether dropping IDEMIX causes the Antarctic summer
sea-ice decline and the rising absorbed solar seen in the TKE-only line after 2190.

Same weather is not expected (the arm starts with a cold atmosphere in 2170), so the comparison is
of decadal means, diff = IDEMIX minus TKE-only, with 2 standard errors from the interannual spread
of each run (no autocorrelation correction), and of the linear trends over 2180-2209.

Quantities as in v350_validation_pairs.py, plus the budget south of 60S: all-sky and clear-sky
absorbed solar at the top of the atmosphere, and the net surface flux.

Usage: idemix_arm_pairs.py     (reval environment; writes data/clim/idemix_arm_pairs.txt)
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
    'idemix_albsn082': dict(root=f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_idemix_albsn082', cache=f'{P}/reval/spinup_cache/idemix_albsn082'),
}
PERIODS = [(2180, 2189), (2190, 2199), (2200, 2209)]
Y0, Y1 = 2171, 2209
ACC, RHO_LF = 3600.0, 3.3355e8
BANDS = [(0, 100), (100, 300), (300, 700), (700, 1500), (1500, 2500), (2500, 6000)]
ZI = np.abs(xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')['nz'].values)
DZ, ZM = np.diff(ZI), 0.5 * (ZI[1:] + ZI[:-1])


def means(root, var, yr):
    """Annual mean of an accumulated field: global and south of 60S."""
    with xr.open_dataset(f'{root}/outdata/oifs/atm_remapped_1m_{var}_{yr}-{yr}.nc') as ds:
        x = ds[var].mean([d for d in ds[var].dims if d not in ('lat', 'lon')]).mean('lon')
        w = np.cos(np.deg2rad(ds['lat']))
        s = ds['lat'] < -60
        return float((x * w).sum() / w.sum()), float((x * w).where(s).sum() / w.where(s).sum())


def table(tag):
    r = RUNS[tag]
    s = pd.read_csv(f'{REPO}/data/clim/series_{tag}.csv').set_index('year')
    rows = []
    for y in range(Y0, Y1 + 1):
        f = {v: means(r['root'], v, y) for v in ('tsr', 'ttr', 'tsrc', 'ssr', 'str', 'sshf', 'slhf', 'sf')}
        g = {v: f[v][0] / ACC for v in f}
        so = {v: f[v][1] / ACC for v in f}
        prof = np.load(f"{r['cache']}/hovm_temp_{y}.npy")
        am = np.load(f"{r['cache']}/amoc_{y}.npy")
        row = dict(year=y, toa=g['tsr'] + g['ttr'], sfc=g['ssr'] + g['str'] + g['sshf'] + g['slhf'] - RHO_LF * g['sf'],
                   asr=g['tsr'], olr=-g['ttr'], asr_so=so['tsr'], asrc_so=so['tsrc'],
                   sfc_so=so['ssr'] + so['str'] + so['sshf'] + so['slhf'] - RHO_LF * so['sf'],
                   amoc26=am[0], amoc_max=am[1])
        n = len(prof)
        for a, b in BANDS:
            m = (ZM[:n] >= a) & (ZM[:n] < b)
            row[f'T{a}-{b}'] = float(np.sum(prof[m] * DZ[:n][m]) / np.sum(DZ[:n][m]))
        rows.append(row)
    d = pd.DataFrame(rows).set_index('year')
    return d.join(s.drop(columns=['exp', 'tsr_raw', 'ttr_raw'], errors='ignore'))


ROWS = [('toa', 'net TOA', 'W/m2', 1), ('sfc', 'net surface', 'W/m2', 1), ('asr', 'absorbed solar', 'W/m2', 1),
        ('olr', 'outgoing longwave', 'W/m2', 1),
        ('asr_so', 'absorbed solar 60-90S', 'W/m2', 1), ('asrc_so', '  clear-sky 60-90S', 'W/m2', 1),
        ('sfc_so', 'net surface 60-90S', 'W/m2', 1),
        ('t2m_glob', 'T2m global', 'C', 1), ('t2m_nh', 'T2m NH', 'C', 1), ('t2m_sh', 'T2m SH', 'C', 1),
        ('t2m_6090n', 'T2m 60-90N', 'C', 1), ('t2m_6090s', 'T2m 60-90S', 'C', 1),
        ('t2m_land', 'T2m land', 'C', 1), ('t2m_ocean', 'T2m ocean', 'C', 1),
        ('sivoln_m04', 'NH ice volume April', '1e3 km3', 1e-3), ('sivoln_m09', 'NH ice volume September', '1e3 km3', 1e-3),
        ('sivols_m09', 'SH ice volume September', '1e3 km3', 1e-3), ('sivols_m02', 'SH ice volume February', '1e3 km3', 1e-3),
        ('siextentn_m03', 'NH ice extent March', '1e6 km2', 1), ('siextentn_m09', 'NH ice extent September', '1e6 km2', 1),
        ('siextents_m09', 'SH ice extent September', '1e6 km2', 1), ('siextents_m02', 'SH ice extent February', '1e6 km2', 1),
        ('amoc26', 'AMOC 26.5N', 'Sv', 1), ('amoc_max', 'AMOC max 20-60N', 'Sv', 1)] + \
       [(f'T{a}-{b}', f'T {a}-{b} m', 'C', 1) for a, b in BANDS]


def trend(y, x):
    """Slope per decade and 2 standard errors."""
    ok = np.isfinite(x)
    y, x = y[ok], x[ok]
    A = np.vstack([y - y.mean(), np.ones_like(y)]).T
    c, res = np.linalg.lstsq(A, x, rcond=None)[:2]
    se = np.sqrt(res[0] / (len(y) - 2) / np.sum((y - y.mean()) ** 2)) if len(res) else np.nan
    return 10 * c[0], 20 * se


A, B = table('tke_albsn082'), table('idemix_albsn082')
out = ['PICAL_crunveg_idemix_albsn082 (TKE + IDEMIX) against PICAL_crunveg_tke_albsn082 (TKE only), branch point end of 2169, albsn 0.82 in both.',
       'Decadal means of each run, diff = IDEMIX minus TKE-only, +- 2 s.e. of the difference; * marks |diff| > 2 s.e.', '']
for p0, p1 in PERIODS:
    out.append(f'{p0}-{p1}')
    out.append(f"{'quantity':26s} {'unit':8s} {'TKE':>10s} {'IDEMIX':>10s} {'diff':>9s} {'+- 2 s.e.':>10s}")
    for key, name, unit, sc in ROWS:
        a, b = A.loc[p0:p1, key].values * sc, B.loc[p0:p1, key].values * sc
        d = np.nanmean(b) - np.nanmean(a)
        se = 2 * np.sqrt(np.nanvar(a, ddof=1) / len(a) + np.nanvar(b, ddof=1) / len(b))
        out.append(f"{name:26s} {unit:8s} {np.nanmean(a):10.3f} {np.nanmean(b):10.3f} {d:+9.3f} {se:10.3f} {'*' if abs(d) > se else ''}")
    out.append('')
out.append('Linear trends over 2180-2209, per decade, +- 2 s.e.; * marks a trend larger than its 2 s.e.')
out.append(f"{'quantity':26s} {'unit':8s} {'TKE':>20s} {'IDEMIX':>20s}")
yy = np.arange(2180, 2210, dtype=float)
for key, name, unit, sc in ROWS:
    ta, tb = trend(yy, A.loc[2180:2209, key].values * sc), trend(yy, B.loc[2180:2209, key].values * sc)
    out.append(f"{name:26s} {unit:8s} {ta[0]:+9.3f} +- {ta[1]:5.3f} {'*' if abs(ta[0]) > ta[1] else ' '} "
               f"{tb[0]:+9.3f} +- {tb[1]:5.3f} {'*' if abs(tb[0]) > tb[1] else ' '}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/idemix_arm_pairs.txt', 'w').write(txt + '\n')
pd.concat([A.assign(run='tke_albsn082'), B.assign(run='idemix_albsn082')]).to_csv(f'{REPO}/data/clim/idemix_arm_pairs_annual.csv', float_format='%.6g')
