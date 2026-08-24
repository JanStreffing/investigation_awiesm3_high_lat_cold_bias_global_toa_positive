"""S4 REMOVAL, coupled, on the current stack: 11V against 11R (1990) and 11W against 11Q (1850).

THE TEST.  S4 (``RCL_INPPMIN`` 70000->50000) was adopted in round 28 on an AMIP Siberian
JJA gain of +1.194* K.  Three charges were laid against it on 2026-08-23/24: coupled, the
11G-11E pair reverses the Siberian sign (-0.365* K); its Southern Ocean benefit is purely
equatorward (-2.34* at 60-45S against +0.77 ns at 90-60S, ocean-masked); and that
equatorward gain is the largest single contributor to the sub-Antarctic overcorrection
that now blocks every SH-wide brightening lever.

11G-11E is four rounds old and sits on a colder, darker base.  These two arms price the
REMOVAL where it would actually be applied, on top of LX4 and RSBLB.  Superposition in
this campaign has been wrong in sign twice, so it is measured rather than inferred.

PRE-REGISTERED FALSIFIERS, from the runscript headers, fixed before the runs finished:

  KEEP S4 if Siberian JJA goes BACKWARDS on removal, i.e. 11V-11R negative.  That would
  mean the 11G-11E reversal does not transfer to this base.
  RETIRE S4 if Siberian JJA is flat or improving, sub-Antarctic DJF SW CRE moves back
  toward zero from -3.70, and the polar band is unchanged.

  Watch item: net TOA.  S4 cools, so removing it warms.  That helps the cold bias and
  spends the energy budget, and 11R sits at -0.027 against its control.

WHY THE 1850 TWIN EXISTS.  The energy question is a PI question and a 1990 arm cannot
answer it: 11R sits at +0.547 over 1352-59 where 11Q sits at -0.079.  The PI controls do
NOT carry a standing +0.88 W/m2 offset, a number this campaign quoted for a week.  They
start near zero and DRIFT, 11N -0.057 -> +0.941 and 11Q -0.079 -> +0.998 over 40 years.
The drift, not the level, is the thing to read, and it is driven by SH September ice
falling ~2e6 km2 in every arm.  So dTOA/dt and d(Sep ice)/dt are scored here explicitly.

TRAPS THIS AVOIDS.
  IFS accumulated fluxes are J/m2 over the output step; divide by 3600 or every radiative
  number is ~3600x too large.
  Band means over the Southern Ocean MUST be ocean-masked.  A 90-60S band average folds
  in the Antarctic ice sheet and halved the measured polar signal once already.
  Sea-ice AREA is not EXTENT.  Both are reported, separately labelled.
  Thresholds come from the CONTROL's own interannual scatter, 1.96*sd*sqrt(2/n), the
  paired form.  A spin-up is drifting rather than equilibrated, so that scatter is an
  upper bound on what two means can resolve, which errs conservative.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
ACC = 3600.0
FLUX = ('tsr', 'ttr', 'tsrc', 'ttrc', 'ssr', 'str')

PAIRS = [('1990', ('11R +S4', f'{R092}/11R'), ('11V -S4', f'{R092}/11V')),
         ('1850', ('11Q +S4', f'{R092}/11Q'), ('11W -S4', f'{R092}/11W'))]

SIB = (55.0, 75.0, 60.0, 180.0)
DJF, JJA = [12, 1, 2], [6, 7, 8]
SEP = 8


def years_avail(root, var):
    D = f'{root}/outdata/oifs'
    if not os.path.isdir(D):
        return []
    ys = []
    for f in os.listdir(D):
        if f.startswith(f'atm_remapped_1m_{var}_') and f.endswith('.nc'):
            try:
                ys.append(int(f.split('_')[-1].split('-')[0]))
            except ValueError:
                pass
    return sorted(ys)


def load_year(root, var, y):
    D = f'{root}/outdata/oifs'
    div = ACC if var in FLUX else 1.0
    f = f'{D}/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(f):
        return None, None, None
    with xr.open_dataset(f, decode_times=False) as d:
        name = var if var in d.data_vars else list(d.data_vars)[-1]
        a = d[name].values / div
        lat, lon = d['lat'].values, d['lon'].values
    return (a, lat, lon) if a.shape[0] == 12 else (None, None, None)


def load_t925(root, y):
    f = f'{root}/outdata/oifs/atm_remapped_1m_pl_t_{y}-{y}.nc'
    if not os.path.exists(f):
        return None
    with xr.open_dataset(f, decode_times=False) as d:
        lv = d['pressure_levels'].values
        k = int(np.argmin(np.abs(lv - 92500.0)))
        return d['t'].values[:, k, :, :]


def fesom_sep(root, var, y):
    g = glob.glob(f'{root}/outdata/fesom/{var}.fesom.{y}.nc')
    if not g:
        return np.nan
    with xr.open_dataset(g[0], decode_times=False) as d:
        return float(d[var].values.ravel()[SEP])


def gmean(f2d, lat):
    return float(np.average(f2d.mean(axis=1), weights=np.cos(np.deg2rad(lat))))


def boxmean(f2d, lat, lon, box, mask=None):
    la0, la1, lo0, lo1 = box
    ys = (lat >= la0) & (lat <= la1)
    xs = ((lon % 360) >= lo0) & ((lon % 360) <= lo1)
    sub = f2d[np.ix_(ys, xs)]
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape).copy()
    if mask is not None:
        w = np.where(mask[np.ix_(ys, xs)], w, 0.0)
    return float(np.average(sub, weights=w)) if w.sum() else np.nan


def zband(f2d, lat, a, b, mask=None):
    sel = (lat >= a) & (lat < b)
    sub = f2d[sel, :]
    w = np.broadcast_to(np.cos(np.deg2rad(lat[sel]))[:, None], sub.shape).copy()
    if mask is not None:
        w = np.where(mask[sel, :], w, 0.0)
    return float(np.average(sub, weights=w)) if w.sum() else np.nan


def metrics_for_year(root, y, lsm, lat, lon):
    t2m, _, _ = load_year(root, '2t', y)
    tsr, _, _ = load_year(root, 'tsr', y)
    ttr, _, _ = load_year(root, 'ttr', y)
    tsrc, _, _ = load_year(root, 'tsrc', y)
    if t2m is None or tsr is None or ttr is None or tsrc is None:
        return None
    land, ocean = lsm >= 0.5, lsm < 0.5
    d = [m - 1 for m in DJF]
    j = [m - 1 for m in JJA]
    cre = (tsr - tsrc)
    net = (tsr + ttr)
    o = {
        'net TOA [W/m2]': gmean(net.mean(axis=0), lat),
        'tropics net TOA [W/m2]': zband(net.mean(axis=0), lat, -30, 30),
        'global T2m [C]': gmean(t2m.mean(axis=0), lat) - 273.15,
        'all-land ANN T2m [C]': zband(t2m.mean(axis=0), lat, -90, 90, land) - 273.15,
        'Arctic land JJA T2m [C]': zband(t2m[j].mean(axis=0), lat, 60, 90, land) - 273.15,
        'Siberia JJA T2m [C]': boxmean(t2m[j].mean(axis=0), lat, lon, SIB, land) - 273.15,
        'Siberia DJF T2m [C]': boxmean(t2m[d].mean(axis=0), lat, lon, SIB, land) - 273.15,
        'polar SW CRE DJF [W/m2]': zband(cre[d].mean(axis=0), lat, -90, -60, ocean),
        'subAnt SW CRE DJF [W/m2]': zband(cre[d].mean(axis=0), lat, -60, -45, ocean),
        'SO SW CRE ANN [W/m2]': zband(cre.mean(axis=0), lat, -65, -45, ocean),
    }
    t925 = load_t925(root, y)
    if t925 is not None:
        inv = (t925 - t2m)[d].mean(axis=0)
        o['Arctic DJF inversion [K]'] = zband(inv, lat, 60, 90, land)
    sst, _, _ = load_year(root, 'sst', y)
    if sst is not None:
        o['SO SST [C]'] = zband(sst.mean(axis=0), lat, -65, -45, ocean) - 273.15
    o['SH Sep ice AREA [1e6 km2]'] = fesom_sep(root, 'siareas', y)
    o['SH Sep ice EXTENT [1e6 km2]'] = fesom_sep(root, 'siextents', y)
    return o


def slope(x, v):
    """OLS slope per decade, and its standard error."""
    x = np.asarray(x, float); v = np.asarray(v, float)
    ok = np.isfinite(v)
    x, v = x[ok], v[ok]
    if len(x) < 5:
        return np.nan, np.nan
    b, a = np.polyfit(x, v, 1)
    res = v - (a + b * x)
    se = np.sqrt((res ** 2).sum() / (len(x) - 2) / ((x - x.mean()) ** 2).sum())
    return b * 10.0, se * 10.0


print(__doc__)

for era, (ctag, croot), (ttag, troot) in PAIRS:
    print('=' * 102)
    print(f'\n### {era} pair:  {ttag} minus {ctag}\n')
    ac, at = years_avail(croot, '2t'), years_avail(troot, '2t')
    win = sorted(set(ac) & set(at))
    print(f'  {ctag:9s} {len(ac):3d} yr {ac[0]}-{ac[-1]}     '
          f'{ttag:9s} {len(at):3d} yr {at[0]}-{at[-1]}')
    if not win:
        print('  no overlap yet\n'); continue
    print(f'  matched window: {len(win)} yr ({win[0]}-{win[-1]})'
          + ('' if len(win) >= 30 else '   BELOW the 30-year campaign minimum'))
    ser = {}
    for tag, root in ((ctag, croot), (ttag, troot)):
        lsm, lat, lon = load_year(root, 'lsm', win[0])
        if lsm is not None and lsm.ndim == 3:
            lsm = lsm[0]
        rows, yrs = [], []
        for y in win:
            m = metrics_for_year(root, y, lsm, lat, lon)
            if m:
                rows.append(m); yrs.append(y)
        ser[tag] = (yrs, rows)
        print(f'  {tag:9s} scored {len(rows)} of {len(win)} years')
    keys = list(ser[ctag][1][0].keys())
    cy, crows = ser[ctag]
    ty, trows = ser[ttag]

    print(f'\n  {"metric":28s} {ctag:>10s} {ttag:>10s} {"diff":>11s} {"thr":>9s}')
    verdict = {}
    for k in keys:
        cv = np.array([r.get(k, np.nan) for r in crows], float)
        tv = np.array([r.get(k, np.nan) for r in trows], float)
        if not np.isfinite(cv).any() or not np.isfinite(tv).any():
            continue
        thr = 1.96 * np.nanstd(cv, ddof=1) * np.sqrt(2.0 / np.isfinite(cv).sum())
        a, b = np.nanmean(cv), np.nanmean(tv)
        dd = b - a
        sig = abs(dd) > thr
        verdict[k] = (a, b, dd, sig, thr)
        print(f'  {k:28s} {a:10.3f} {b:10.3f} {dd:+10.3f}{"*" if sig else " "} {thr:9.3f}')

    print(f'\n  DRIFT over the same window, per decade\n')
    print(f'  {"metric":28s} {ctag:>12s} {ttag:>12s} {"d(slope)":>12s}')
    for k in ('net TOA [W/m2]', 'SH Sep ice AREA [1e6 km2]',
              'SH Sep ice EXTENT [1e6 km2]', 'SO SST [C]', 'global T2m [C]'):
        if k not in verdict:
            continue
        bc, sc = slope(cy, [r.get(k, np.nan) for r in crows])
        bt, st = slope(ty, [r.get(k, np.nan) for r in trows])
        if not np.isfinite(bc) or not np.isfinite(bt):
            continue
        dse = np.sqrt(sc ** 2 + st ** 2)
        mark = '*' if abs(bt - bc) > 1.96 * dse else ' '
        print(f'  {k:28s} {bc:+8.3f}+-{sc:.3f} {bt:+8.3f}+-{st:.3f} '
              f'{bt - bc:+10.3f}{mark}')

    print('\n  PRE-REGISTERED FALSIFIER')
    if 'Siberia JJA T2m [C]' in verdict:
        a, b, dd, sig, thr = verdict['Siberia JJA T2m [C]']
        if dd < 0 and sig:
            print(f'    Siberia JJA {dd:+.3f}* BACKWARDS and resolved -> KEEPS S4')
        elif dd < 0:
            print(f'    Siberia JJA {dd:+.3f} backwards but within +-{thr:.3f} -> not a keep')
        else:
            print(f'    Siberia JJA {dd:+.3f}{"*" if sig else ""} flat or improving '
                  f'-> the keep-S4 falsifier did NOT fire')
    for k in ('subAnt SW CRE DJF [W/m2]', 'polar SW CRE DJF [W/m2]', 'net TOA [W/m2]'):
        if k in verdict:
            a, b, dd, sig, thr = verdict[k]
            print(f'    {k:28s} {a:+8.3f} -> {b:+8.3f}  ({dd:+.3f}'
                  f'{", resolved" if sig else ", ns"})')
    print()
