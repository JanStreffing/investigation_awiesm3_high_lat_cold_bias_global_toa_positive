"""The three-arm land-surface test: 11G / 11I / 11J, scored as two pairs.

THE DESIGN, as pre-registered in the runscript headers:

    11G  old slt, no Raupach     the adopted stack (S4 = RCL_INPPMIN 50k)   50 yr
    11I  _v2 slt, no Raupach     the soil-map repair, priced on its own     50 yr
    11J  _v2 slt, Raupach on     roughness, priced on top of the repair     30 yr

    11I - 11G  = the land-surface repair
    11J - 11I  = Raupach canopy roughness

PRE-REGISTERED READING FOR 11I - 11G (thresholds +-0.37 K Siberian DJF soil,
+-0.61 K DJF T2m at 30 coupled years):
    nothing resolves   -> the repair is free, adopt and move on
    Siberia warms      -> the old map was part of the cold bias
    Siberia cools      -> the repair costs land temperature and the land campaign's
                          numbers need restating against it

PRE-REGISTERED READING FOR 11J - 11I.  Raupach is a WINTER lever: FracHVeg reports a
summergreen boreal canopy as nearly absent in winter (Siberian cvh 0.239 Aug -> 0.091
Jan in 11E), handing ~91 % of the winter grid box the bare-soil 0.013 m in polar night,
whereas Raupach makes a leafless canopy the roughest state.  Offline the grid-box
roughness goes DJF 0.058 -> 0.140 m and JJA 0.190 -> 0.120 m.  So:
    DJF soil warms, JJA inside noise -> adopt
    DJF soil does not resolve        -> the roughness route is weaker than the campaign
                                        regression (+0.4 cvh -> +1.2 K) implied
    JJA cools beyond its threshold   -> the summer cost is real and must be priced
NOT a forest fix: BNS clears every climate gate already and is limited by C3 grass
competition, which roughness does not touch.  Do not read tree cover as the verdict.

HONESTY ABOUT WHAT EACH PAIR ISOLATES.

  11J - 11I IS a one-variable pair.  Both arms load the SAME binary (md5 ff9b5507) and
  differ only by `ifraupachz0` in guess.ins plus the matching ECE_CPL_LPJG_Z0 in fort.4.
  The Raupach code is compiled into both and gated at run time precisely so that this
  comparison has no binary difference in it.

  11I - 11G IS NOT.  The header calls it "ONE input file", which was true of the
  original design but is NOT true of what ran.  Three things differ:
      1. lpjg_slt_suffix _v2      -- 302 simulated cells moved off soil code 6
      2. the LPJ-GUESS gridlist   -- L096 -> TCO95-land (mask fix 168a98b)
      3. the binary               -- 11G md5 8c5ab467, 11I ff9b5507, which adds the
                                     cold-start path; 307 cells that have no spin-up
                                     state cold-start in 11I instead of aborting
  2 and 3 are not separable: the new gridlist is what ADDS those cells, and the
  cold-start fix is what lets them run.  So read 11I - 11G as "the land-surface repair
  bundle", not as the soil map alone.  A clean one-variable soil-map test would need an
  11I-prime on the L096 gridlist, which has not been run.

WINDOWS.  1370-79 is excluded from 11I - 11G: the repaired CORE3 soil map slipped into
11G's leg 3 alone (organic cell count 233 there against 547 in its other four legs), so
in that decade 11G is partly _v2 itself.  That leaves 1350-69 + 1380-99 = 40 years.
11J - 11I uses 1350-79, all 30 years 11J produced before it lost leg 4 to a NetCDF/HDF
write failure at finalisation.

THRESHOLDS FIRST, from the CONTROL arm's own interannual scatter over the same window,
1.96*sd*sqrt(2/n).  A coupled spin-up is drifting rather than equilibrated, so its
scatter is an upper bound on what a pair of means can resolve -- the conservative
direction.

TRAP: IFS TOA and surface fluxes are accumulated J/m^2 over the output step and must be
divided by the accumulation period.  load_year does it for the flux variables only.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
ACC = 3600.0

ARM = {
    '11G': f'{R092}/Tuning_test_11G_inppmin50k',
    '11I': f'{R092}/Tuning_test_11I_v2soil',
    '11J': f'{R092}/Tuning_test_11J_v2soil_raupach',
}

# (control, treatment, window, label)
PAIRS = [
    ('11G', '11I', list(range(1350, 1370)) + list(range(1380, 1400)),
     'the land-surface repair bundle (_v2 slt + TCO95-land gridlist + cold-start)'),
    ('11I', '11J', list(range(1350, 1380)),
     'Raupach canopy roughness, one namelist flag, one shared binary'),
]

SIB = (55.0, 75.0, 60.0, 180.0)
SO = (-65.0, -45.0)
DJF, JJA = [12, 1, 2], [6, 7, 8]
FLUX = ('tsr', 'ttr', 'tsrc', 'ttrc', 'ssr', 'str')


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
    for f in (f'{D}/atm_remapped_1m_{var}_{y}-{y}.nc',
              f'{D}/atm_remapped_1m_{var}_1m_{y}-{y}.nc'):
        if os.path.exists(f):
            with xr.open_dataset(f, decode_times=False) as d:
                a = d[var].values / div
                lat, lon = d['lat'].values, d['lon'].values
            return (a, lat, lon) if a.shape[0] == 12 else (None, None, None)
    return None, None, None


def gmean(f2d, lat):
    return float(np.average(f2d.mean(axis=1), weights=np.cos(np.deg2rad(lat))))


def boxmean(f2d, lat, lon, box, lsm=None):
    la0, la1, lo0, lo1 = box
    ys = (lat >= la0) & (lat <= la1)
    xs = ((lon % 360) >= lo0) & ((lon % 360) <= lo1)
    sub = f2d[np.ix_(ys, xs)]
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape).copy()
    if lsm is not None:
        w = np.where(lsm[np.ix_(ys, xs)] >= 0.5, w, 0.0)
    return float(np.average(sub, weights=w)) if w.sum() else np.nan


def zband(f2d, lat, a, b):
    sel = (lat >= a) & (lat < b)
    return float(np.average(f2d[sel, :].mean(axis=1),
                            weights=np.cos(np.deg2rad(lat[sel]))))


def zband_ocean(f2d, lat, a, b, lsm):
    sel = (lat >= a) & (lat < b)
    w = np.broadcast_to(np.cos(np.deg2rad(lat[sel]))[:, None], f2d[sel, :].shape).copy()
    w = np.where(lsm[sel, :] < 0.5, w, 0.0)
    return float(np.average(f2d[sel, :], weights=w)) if w.sum() else np.nan


def metrics_for_year(root, y, lsm, lat, lon):
    t2m, _, _ = load_year(root, '2t', y)
    tsr, _, _ = load_year(root, 'tsr', y)
    ttr, _, _ = load_year(root, 'ttr', y)
    if t2m is None or tsr is None or ttr is None:
        return None
    tsrc, _, _ = load_year(root, 'tsrc', y)
    stl2, _, _ = load_year(root, 'stl2', y)
    d = [m - 1 for m in DJF]
    j = [m - 1 for m in JJA]
    out = {
        'net TOA [W/m2]': gmean((tsr + ttr).mean(axis=0), lat),
        'SO SW CRE [W/m2]': (zband((tsr - tsrc).mean(axis=0), lat, *SO)
                             if tsrc is not None else np.nan),
        'SO net TOA [W/m2]': zband((tsr + ttr).mean(axis=0), lat, *SO),
        'global T2m [C]': gmean(t2m.mean(axis=0), lat) - 273.15,
        'Siberia JJA T2m [C]': boxmean(t2m[j].mean(axis=0), lat, lon, SIB, lsm) - 273.15,
        'Siberia DJF T2m [C]': boxmean(t2m[d].mean(axis=0), lat, lon, SIB, lsm) - 273.15,
    }
    sst, _, _ = load_year(root, 'sst', y)
    ci, _, _ = load_year(root, 'ci', y)
    if tsrc is not None:
        out['SO clear-sky SW [W/m2]'] = zband(tsrc.mean(axis=0), lat, *SO)
    if sst is not None:
        out['SO SST [C]'] = zband_ocean(sst.mean(axis=0), lat, *SO, lsm) - 273.15
    if ci is not None:
        out['SO sea ice [%]'] = zband_ocean(ci.mean(axis=0), lat, *SO, lsm) * 100.0
    if stl2 is not None:
        out['Siberia DJF soil [C]'] = boxmean(stl2[d].mean(axis=0), lat, lon, SIB, lsm) - 273.15
        out['Siberia JJA soil [C]'] = boxmean(stl2[j].mean(axis=0), lat, lon, SIB, lsm) - 273.15
    return out


def series_for(tag, years):
    root = ARM[tag]
    lsm, lat, lon = load_year(root, 'lsm', years[0])
    if lsm is not None and lsm.ndim == 3:
        lsm = lsm[0]
    rows = []
    for y in years:
        m = metrics_for_year(root, y, lsm, lat, lon)
        if m:
            rows.append(m)
    return rows


print(__doc__)
print('=' * 100)
print('\nYEARS ON DISK\n')
avail = {}
for tag, root in ARM.items():
    ys = years_avail(root, '2t')
    avail[tag] = ys
    print(f'  {tag}  {len(ys):3d} years' + (f'  {ys[0]}-{ys[-1]}' if ys else '  NO OUTPUT'))

ceres_so = None
try:
    with xr.open_dataset(CERESF) as cds:
        clat = cds['lat'].values
        csel = (clat >= SO[0]) & (clat < SO[1])
        ceres_so = float(np.average(
            cds['toa_cre_sw_clim'].values.mean(axis=0)[csel, :].mean(axis=1),
            weights=np.cos(np.deg2rad(clat[csel]))))
except Exception as exc:
    print(f'\n  (CERES anchor unavailable: {exc})')

for ctl, trt, window, label in PAIRS:
    usable = [y for y in window if y in avail[ctl] and y in avail[trt]]
    print('\n' + '=' * 100)
    print(f'\n{trt} MINUS {ctl}  --  {label}')
    if not usable:
        print('  no overlapping years; skipped')
        continue
    print(f'  window: {len(usable)} years, {usable[0]}-{usable[-1]}'
          + ('' if len(usable) >= 30 else '   WARNING: below the 30-year minimum'))
    missing = [y for y in window if y not in usable]
    if missing:
        print(f'  requested but unavailable: {missing[0]}-{missing[-1]} '
              f'({len(missing)} years)')

    sc, st = series_for(ctl, usable), series_for(trt, usable)
    print(f'  scored {len(sc)} / {len(st)} years ({ctl} / {trt})')
    if not sc or not st:
        continue
    keys = [k for k in sc[0] if k in st[0]]

    print(f'\n  detection thresholds from {ctl} interannual scatter, '
          f'1.96*sd*sqrt(2/n):\n')
    thr = {}
    for k in keys:
        v = np.array([r[k] for r in sc])
        thr[k] = 1.96 * v.std(ddof=1) * np.sqrt(2.0 / len(v))

    print(f'  {"metric":24s} {ctl:>10s} {trt:>10s} {"diff":>11s} {"thr":>9s}')
    verdict = {}
    for k in keys:
        a = float(np.mean([r[k] for r in sc]))
        b = float(np.mean([r[k] for r in st]))
        d = b - a
        sig = abs(d) > thr[k]
        verdict[k] = (a, b, d, sig)
        print(f'  {k:24s} {a:10.3f} {b:10.3f} {d:+10.3f}{"*" if sig else " "} {thr[k]:9.3f}')

    if ceres_so is not None and 'SO SW CRE [W/m2]' in verdict:
        a, b, _, _ = verdict['SO SW CRE [W/m2]']
        print(f'\n  CERES SO SW CRE {ceres_so:.2f}:  {ctl} is {a - ceres_so:+.2f} from it, '
              f'{trt} {b - ceres_so:+.2f}')

    print('\n  PRE-REGISTERED READING')
    res = [k for k in keys if verdict[k][3]]
    if not res:
        print('    nothing resolves against its own threshold')
    else:
        for k in res:
            print(f'    RESOLVED  {k:24s} {verdict[k][2]:+.3f}  (thr {thr[k]:.3f})')
    for k in ('Siberia DJF soil [C]', 'Siberia JJA soil [C]',
              'Siberia DJF T2m [C]', 'Siberia JJA T2m [C]', 'net TOA [W/m2]'):
        if k in verdict:
            a, b, d, s = verdict[k]
            print(f'    {k:24s} {d:+8.3f}  {"resolved" if s else "ns"}')
print()
