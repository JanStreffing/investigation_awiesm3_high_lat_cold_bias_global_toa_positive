"""Pack-ice drift, air-ice stress and Southern Ocean deep convection, per run and period.

For runs branched from the PICAL_crunveg_tke_albsn082 restart of 2169. Per year:
  SH pack, JJA (south of 55S, daily ice concentration >= 0.8): mean drift speed and mean northward
    drift from the daily uice, vice; mean magnitude of the monthly air-ice and ice-ocean stress
    (monthly concentration >= 0.8).
  NH pack, JFM (north of 70N): the same, without the northward component.
  Area of the open ocean south of 55S whose September mixed layer (MLD2) is deeper than 1000 m.
Period means, and the standard deviation of the annual values in brackets for the first run listed.

Usage: ice_speed_and_so_convection.py runs periods outname
  e.g. z0def,z0floor,z0form 2171-2179,2180-2189 z0_ice_speed   (reval environment; writes data/clim/<outname>.txt)
"""
import os, sys
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
EXP = lambda k: 'PICAL_crunveg_tke_albsn082' if k == 'control' else f'PICAL_crunveg_{k}'
SEL = sys.argv[1].split(',')
PER = [tuple(int(v) for v in p.split('-')) for p in sys.argv[2].split(',')]
OUT = sys.argv[3]
md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lat = md['lat'].values; area = md['nod_area'].values[0]
od = (md['ulevels_nod2D'].values == 1) & (lat < -55)
COLS = ['SH speed', 'SH vnorth', 'SH tau_ai', 'SH tau_io', 'NH speed', 'NH tau_ai', 'NH tau_io', 'conv area']


def year(k, y):
    d = f'{P}/runtime/awiesm3-v3.4/{EXP(k)}/outdata/fesom'
    o = lambda v: xr.open_dataset(f'{d}/{v}.fesom.{y}.nc')[v]
    a, u, v = o('a_ice'), o('uice'), o('vice')
    mon = a['time'].dt.month.values
    am = a.groupby('time.month').mean('time').values                      # (12, nod)
    r = {}
    for hem, reg, months in (('SH', lat < -55, (6, 7, 8)), ('NH', lat > 70, (1, 2, 3))):
        t = np.isin(mon, months); idx = np.where(reg)[0]
        aa, uu, vv = a.values[t][:, idx], u.values[t][:, idx], v.values[t][:, idx]
        w = np.where((aa >= 0.8) & np.isfinite(uu) & np.isfinite(vv), area[idx][None, :], 0.0)
        r[f'{hem} speed'] = 100 * np.nansum(np.hypot(uu, vv) * w) / w.sum()
        if hem == 'SH':
            r['SH vnorth'] = 100 * np.nansum(vv * w) / w.sum()
        for nm, vx, vy in (('tau_ai', 'atmice_x', 'atmice_y'), ('tau_io', 'iceoce_x', 'iceoce_y')):
            tx, ty = o(vx).values, o(vy).values
            num = den = 0.0
            for m in months:
                tm = np.hypot(tx[m - 1][idx], ty[m - 1][idx])
                wm = np.where((am[m - 1][idx] >= 0.8) & np.isfinite(tm), area[idx], 0.0)
                num += np.nansum(tm * wm); den += wm.sum()
            r[f'{hem} {nm}'] = 1e3 * num / den
    m = np.abs(o('MLD2').isel(time=8).values)
    r['conv area'] = np.sum(area[od & (m > 1000)]) / 1e12
    return [r[c] for c in COLS]


out = ['Pack ice (concentration >= 0.8): drift speed and northward drift [cm/s], air-ice and ice-ocean stress magnitude [mN/m2];',
       'SH south of 55S in JJA, NH north of 70N in JFM. conv area: open ocean south of 55S with September MLD2 > 1000 m [1e6 km2].',
       f'In brackets: standard deviation of the annual values of {SEL[0]}.', '']
for p0, p1 in PER:
    out += [f'{p0}-{p1}', f"{'run':10s}" + ''.join(f'{c:>12s}' for c in COLS)]
    for i, k in enumerate(SEL):
        ys = [y for y in range(p0, p1 + 1) if os.path.isfile(f'{P}/runtime/awiesm3-v3.4/{EXP(k)}/outdata/fesom/uice.fesom.{y}.nc')]
        if not ys:
            out.append(f'{k:10s}   no output'); continue
        X = np.array([year(k, y) for y in ys])
        out.append(f'{k:10s}' + ''.join(f'{x:12.2f}' for x in X.mean(0)))
        if i == 0:
            out.append(f"{'':10s}" + ''.join(f"{'(' + format(s, '.2f') + ')':>12s}" for s in X.std(0, ddof=1)))
        print(k, p0, p1, flush=True)
    out.append('')
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/{OUT}.txt', 'w').write(txt + '\n')
