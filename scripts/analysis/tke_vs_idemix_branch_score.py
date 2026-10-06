"""Score PICAL_crunveg_gmhemi_tke (cvmix_TKE, IDEMIX off, 2140-2159) against its parent
PICAL_crunveg_gmhemi1800 (cvmix_TKE + IDEMIX, 2120-2139).

The two are SEQUENTIAL, not parallel: the TKE run starts from the parent's 2139-12-31
restart, and no IDEMIX run exists for 2140-2159.  So every number is a 20-year segment of
one trajectory, and a difference between segments mixes the scheme change with whatever
drift the parent was on.  The script therefore prints trends per segment next to segment
means, and the figure shows both segments on one axis.

Surface:     T2m by region, sea-ice volume and extent (series_<run>.csv, climate_annual_series.py)
Subsurface:  global-mean temperature per depth band (spin-up cache, cdo fldmean per level,
             bands averaged by layer thickness), thetaoga, AMOC at 26.5N
Net flux:    TOA (tsr + ttr) and surface (ssr + str + sshf + slhf - rho_w Lf sf), global means
             of the remapped monthly OpenIFS files; ocean heat content tendency from thetaoga

Usage: tke_vs_idemix_branch_score.py   (reval environment; writes data/clim/ and plots/)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
RUNS = {
    'gmhemi1800': dict(root=f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_gmhemi1800', y0=2120, y1=2139,
                       cache=f'{P}/reval/spinup_cache/critical_path', label='TKE + IDEMIX (gmhemi1800)', c='tab:red'),
    'gmhemi_tke': dict(root=f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_gmhemi_tke', y0=2140, y1=2159,
                       cache=f'{P}/reval/spinup_cache/gmhemi_tke', label='TKE only (gmhemi_tke)', c='tab:blue'),
    'tke_albsn082': dict(root=f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082', y0=2160, y1=2209,
                         cache=f'{P}/reval/spinup_cache/tke_albsn082', label='TKE only, albsn 0.82', c='tab:green'),
}
ACC = 3600.0
RHO_LF = 3.3355e8          # rho_water * Lf, sf is metres of water per output step
BANDS = [(0, 100), (100, 300), (300, 700), (700, 1500), (1500, 2500), (2500, 6000)]
A_EARTH = 5.101e14
RHO_CP = 1026.0 * 3990.0

zdiag = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
ZI = np.abs(zdiag['nz'].values)            # 48 interfaces
DZ, ZM = np.diff(ZI), 0.5 * (ZI[1:] + ZI[:-1])


def gmean(root, var, yr):
    with xr.open_dataset(f'{root}/outdata/oifs/atm_remapped_1m_{var}_{yr}-{yr}.nc') as ds:
        x = ds[var].mean('time_counter' if 'time_counter' in ds[var].dims else 'time')
        w = np.cos(np.deg2rad(ds['lat']))
        return float((x * w).sum(('lat', 'lon')) / (w.sum('lat') * ds.sizes['lon']))


def fesom_scalar(root, var, yr):
    fn = f'{root}/outdata/fesom/{var}.fesom.{yr}.nc'
    if not os.path.exists(fn):
        return np.nan
    with xr.open_dataset(fn) as ds:
        return float(np.mean(ds[var].values))


def run_table(tag):
    r = RUNS[tag]
    s = pd.read_csv(f'{REPO}/data/clim/series_{tag}.csv').set_index('year')
    rows = []
    for y in range(r['y0'], r['y1'] + 1):
        f = {v: gmean(r['root'], v, y) / ACC for v in ('tsr', 'ttr', 'ssr', 'str', 'sshf', 'slhf', 'sf')}
        prof = np.load(f"{r['cache']}/hovm_temp_{y}.npy")
        am = np.load(f"{r['cache']}/amoc_{y}.npy")
        row = dict(year=y, toa=f['tsr'] + f['ttr'],
                   sfc=f['ssr'] + f['str'] + f['sshf'] + f['slhf'] - RHO_LF * f['sf'],
                   asr=f['tsr'], olr=-f['ttr'], amoc26=am[0], amoc_max=am[1],
                   thetaoga=fesom_scalar(r['root'], 'thetaoga', y), volo=fesom_scalar(r['root'], 'volo', y))
        n = len(prof)
        for a, b in BANDS:
            m = (ZM[:n] >= a) & (ZM[:n] < b)
            row[f'T{a}-{b}'] = float(np.sum(prof[m] * DZ[:n][m]) / np.sum(DZ[:n][m]))
        rows.append(row)
    d = pd.DataFrame(rows).set_index('year')
    return d.join(s.drop(columns=['exp', 'tsr_raw', 'ttr_raw'], errors='ignore'))


def trend(y, x, per=10.0):
    """Slope per decade and 2 standard errors (no autocorrelation correction)."""
    ok = np.isfinite(x)
    y, x = np.asarray(y)[ok], np.asarray(x)[ok]
    A = np.vstack([y - y.mean(), np.ones_like(y)]).T
    co, res, *_ = np.linalg.lstsq(A, x, rcond=None)
    se = np.sqrt(res[0] / (len(y) - 2) / np.sum((y - y.mean()) ** 2)) if len(res) else np.nan
    return co[0] * per, 2 * se * per


T = {k: run_table(k) for k in RUNS if os.path.exists(f'{REPO}/data/clim/series_{k}.csv')}
allT = pd.concat(T.values())
allT.to_csv(f'{REPO}/data/clim/tke_vs_idemix_branch_annual.csv', float_format='%.6g')

ROWS = [('toa', 'net TOA', 'W/m2', 1), ('sfc', 'net surface (incl. snow enthalpy)', 'W/m2', 1),
        ('asr', 'absorbed solar', 'W/m2', 1), ('olr', 'outgoing longwave', 'W/m2', 1),
        ('t2m_glob', 'T2m global', 'C', 1), ('t2m_nh', 'T2m NH', 'C', 1), ('t2m_sh', 'T2m SH', 'C', 1),
        ('t2m_6090n', 'T2m 60-90N', 'C', 1), ('t2m_6090s', 'T2m 60-90S', 'C', 1),
        ('t2m_land', 'T2m land', 'C', 1), ('t2m_ocean', 'T2m ocean', 'C', 1),
        ('sivoln_m04', 'NH ice volume April', '1e3 km3', 1e-3), ('sivoln_m09', 'NH ice volume September', '1e3 km3', 1e-3),
        ('sivols_m09', 'SH ice volume September', '1e3 km3', 1e-3), ('sivols_m02', 'SH ice volume February', '1e3 km3', 1e-3),
        ('siextentn_m03', 'NH ice extent March', '1e6 km2', 1), ('siextentn_m09', 'NH ice extent September', '1e6 km2', 1),
        ('siextents_m09', 'SH ice extent September', '1e6 km2', 1), ('siextents_m02', 'SH ice extent February', '1e6 km2', 1),
        ('amoc26', 'AMOC 26.5N', 'Sv', 1), ('amoc_max', 'AMOC max 20-60N', 'Sv', 1),
        ('thetaoga', 'ocean mean temperature', 'C', 1)] + \
       [(f'T{a}-{b}', f'T {a}-{b} m', 'C', 1) for a, b in BANDS]

def compare(ka, kb, na, nb, head):
    """Segment kb against its parent ka: mean of the last na / nb years, trend per decade over each
    whole segment, and the step across the branch (first 5 years of kb minus last 5 of ka)."""
    o = list(head)
    o.append(f"{'quantity':36s} {'unit':8s} {'parent mean':>12s} {'branch mean':>12s} {'diff':>8s} {'parent trend/dec':>18s} {'branch trend/dec':>18s} {'step':>8s}")
    a, b = T[ka], T[kb]
    for key, name, unit, sc in ROWS:
        if key not in a:
            continue
        xa, xb = a[key].values * sc, b[key].values * sc
        ma, mb = np.nanmean(xa[-na:]), np.nanmean(xb[-nb:])
        ta, tb = trend(a.index.values, xa), trend(b.index.values, xb)
        step = np.nanmean(xb[:5]) - np.nanmean(xa[-5:])
        o.append(f'{name:36s} {unit:8s} {ma:12.3f} {mb:12.3f} {mb - ma:+8.3f} {ta[0]:+9.3f} ±{ta[1]:6.3f} {tb[0]:+9.3f} ±{tb[1]:6.3f} {step:+8.3f}')
    return o


out = compare('gmhemi1800', 'gmhemi_tke', 15, 15, [
    'PICAL_crunveg_gmhemi1800 (TKE+IDEMIX, 2120-2139) and its continuation PICAL_crunveg_gmhemi_tke (TKE only, 2140-2159).',
    'Sequential segments of one trajectory. Means over the last 15 years of each segment (2125-2139, 2145-2159); trend over',
    'the 20-year segment per decade +- 2 s.e.; step = first 5 years of the branch minus the last 5 years of the parent.'])
if 'tke_albsn082' in T:
    out.append('')
    out += compare('gmhemi_tke', 'tke_albsn082', 10, 15, [
        'PICAL_crunveg_gmhemi_tke (albsn 0.80, 2140-2159) and its continuation PICAL_crunveg_tke_albsn082 (albsn 0.82, 2160-2209).',
        'Parent mean 2150-2159, branch mean 2195-2209; parent trend over 20 years, branch trend over all 50 years;',
        'step = 2160-2164 minus 2155-2159.'])
    d = T['tke_albsn082']
    out.append('')
    out.append('PICAL_crunveg_tke_albsn082 by decade (2160 excluded from the first): means, and for the last line the trend over 2180-2209 per decade +- 2 s.e.')
    keys = [k for k in ('toa', 'sfc', 't2m_glob', 't2m_6090n', 't2m_6090s', 'sivoln_m04', 'sivoln_m09', 'sivols_m09', 'sivols_m02',
                        'siextentn_m09', 'siextents_m02', 'amoc26', 'T0-100', 'T300-700', 'T700-1500') if k in d]
    out.append(f"{'years':12s} " + ' '.join(f'{k:>13s}' for k in keys))
    for a0 in range(2160, 2210, 10):
        s_ = d[(d.index >= max(a0, 2161)) & (d.index < a0 + 10)]
        out.append(f'{max(a0, 2161)}-{a0 + 9}    ' + ' '.join(f'{s_[k].mean():13.4g}' for k in keys))
    l = d[d.index >= 2180]
    out.append(f"{'trend 2180+':12s} " + ' '.join(f'{trend(l.index.values, l[k].values)[0]:+13.3g}' for k in keys))
    out.append(f"{'  +- 2 s.e.':12s} " + ' '.join(f'{trend(l.index.values, l[k].values)[1]:13.3g}' for k in keys))

# ocean heat uptake implied by thetaoga, W/m2 of Earth surface
out.append('')
for k, d in T.items():
    if np.isfinite(d['thetaoga']).sum() > 3:
        tr, se = trend(d.index.values, d['thetaoga'].values, per=1.0)
        vol = np.nanmean(d['volo'])
        f = tr * RHO_CP * vol / (365.25 * 86400) / A_EARTH
        out.append(f"{k:12s} ocean heat uptake from the thetaoga trend: {f:+.3f} ± {se * RHO_CP * vol / (365.25 * 86400) / A_EARTH:.3f} W/m2 (Earth surface); "
                   f"mean net TOA {d['toa'].mean():+.3f}, mean net surface {d['sfc'].mean():+.3f}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/tke_vs_idemix_branch_score.txt', 'w').write(txt + '\n')

# ---------------------------------------------------------------- figure
fig, axs = plt.subplots(3, 3, figsize=(15, 10.5), sharex=True)
panels = [('toa', 'net TOA [W/m2]'), ('t2m_glob', 'T2m global [C]'), ('T700-1500', 'T 700-1500 m [C]'),
          ('sivoln_m04', 'NH ice volume April [km3]'), ('sivoln_m09', 'NH ice volume September [km3]'), ('siextentn_m09', 'NH ice extent September [1e6 km2]'),
          ('sivols_m09', 'SH ice volume September [km3]'), ('sivols_m02', 'SH ice volume February [km3]'), ('siextents_m02', 'SH ice extent February [1e6 km2]')]
for ax, (key, title), lab in zip(axs.flat, panels, 'abcdefghi'):
    for k, d in T.items():
        r = RUNS[k]
        ax.plot(d.index, d[key], color=r['c'], lw=1.0, marker='o', ms=2.5, label=r['label'])
        tr = np.polyfit(d.index, d[key], 1)
        ax.plot(d.index, np.polyval(tr, d.index), color=r['c'], lw=1.6, alpha=0.5)
    for xb in (2139.5, 2159.5)[:len(T) - 1]:
        ax.axvline(xb, color='0.3', ls=':', lw=0.8)
    if key in ('toa', 'sfc'):
        ax.axhline(0, color='0.6', lw=0.6)
    ax.set_title(f'({lab}) {title}', loc='left', fontsize=11)
    ax.grid(color='0.9', lw=0.5)
axs[0, 0].legend(fontsize=9, loc='best')
for ax in axs[-1]:
    ax.set_xlabel('year')
fig.suptitle('IDEMIX off at 2140, albsn 0.80 to 0.82 at 2160: annual means and linear fits per segment', fontweight='bold')
fig.tight_layout()
fig.savefig(f'{REPO}/plots/tke_vs_idemix_branch_score.png', dpi=160)
