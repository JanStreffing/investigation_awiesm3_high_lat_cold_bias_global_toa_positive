"""S4 coupled: 11G against 11E, on the years where it is a one-variable pair.

THE TEST.  Round 28 scored S4 (``RCL_INPPMIN`` 70000->50000) the best AMIP arm on all
three constraints at once: it is the only one that improves the Southern Ocean surface
while landing global net TOA near the piControl target, at no cost in either Siberian
season or the tropics, and it is one namelist number.  Every coupled arm ever run --
10A, 10B, 11A-11F -- uses 70000.  11G is 11E with that single line changed, against a
run already complete at 50 years, so it needs no new control.

The AMIP number should be a LOWER bound: round 27 measured the coupled amplification at
1.40x, and it comes from sea-ice albedo rather than cloud, so a lever that cools the SO
*surface* is exactly the kind that should amplify.

WHY 1370-79 IS EXCLUDED, and why this is not optional.  Measured from the staged
``slt_TCO95.nc`` in each leg's work directory, organic (type 6) cell counts:

    years      11E    11G
    1350-59    547    547     matched
    1360-69    547    547     matched
    1370-79    547    233     <-- the repaired CORE3 soil map slipped into 11G alone
    1380-89    547    547     matched
    1390-99    547    547     matched (11G leg in flight)

For that one decade 11G differs from 11E in TWO things: RCL_INPPMIN and 314 cells moved
off soil code 6, which soil.cpp treats on its own branch.  Using "the last 30 years" --
the campaign's usual coupled window, and what coupled_r11_eval.py would do -- lands
exactly on 1370-99 and silently includes it.  This script uses the matched set instead.

All five legs of both runs loaded the same LPJ-GUESS binary (md5 8c5ab467), so there is
no binary discontinuity alongside the soil one.

THRESHOLDS FIRST.  Every difference is scored against a detection threshold computed
from 11E's own interannual scatter over the same window, 1.96*sd*sqrt(2/n).  A coupled
spin-up is drifting, not equilibrated, so the scatter is an upper bound on what a pair
of means can resolve -- which is the conservative direction.

A TRAP THIS AVOIDS.  IFS TOA fluxes are accumulated J/m^2 over the output step and must
be divided by the accumulation period; forgetting it makes every radiative number ~3600x
too large.  ``load_year`` below does it for the flux variables only.
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

PAIR = [('11E base', f'{R092}/Tuning_test_11E_swemin15_K1'),
        ('11G +S4', f'{R092}/Tuning_test_11G_inppmin50k')]

# The matched decades. 1370-79 is deliberately absent; see the docstring.
CLEAN = list(range(1350, 1370)) + list(range(1380, 1390))
DIRTY = list(range(1370, 1380))

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
    """(12, nlat, nlon) for one year, de-accumulated, plus lat/lon."""
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
    w = np.cos(np.deg2rad(lat))
    return float(np.average(f2d.mean(axis=1), weights=w))


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
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(f2d[sel, :].mean(axis=1), weights=w))


def zband_ocean(f2d, lat, a, b, lsm):
    """Zonal band mean over OCEAN points only.

    The surface terms have to be ocean-masked: a band mean of SST or sea-ice
    cover that includes the Antarctic Peninsula and the southern tips of the
    continents is not a Southern Ocean number.
    """
    sel = (lat >= a) & (lat < b)
    w = np.broadcast_to(np.cos(np.deg2rad(lat[sel]))[:, None], f2d[sel, :].shape).copy()
    w = np.where(lsm[sel, :] < 0.5, w, 0.0)
    return float(np.average(f2d[sel, :], weights=w)) if w.sum() else np.nan


def metrics_for_year(root, y, lsm, lat, lon):
    """The scored quantities for a single year, or None if output is incomplete."""
    t2m, _, _ = load_year(root, '2t', y)
    if t2m is None:
        return None
    tsr, _, _ = load_year(root, 'tsr', y)
    ttr, _, _ = load_year(root, 'ttr', y)
    tsrc, _, _ = load_year(root, 'tsrc', y)
    stl2, _, _ = load_year(root, 'stl2', y)
    if tsr is None or ttr is None:
        return None
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
    # The terms the coupled amplification actually lives in. Round 27 measured
    # 1.40x for DMS and found the coupled CRE response SMALLER than AMIP's --
    # the cloud adjusts away once the surface cools, and the amplification comes
    # from sea-ice albedo, which is CLEAR-SKY and invisible to CRE. So scoring an
    # SO lever on CRE alone measures the part that does not survive coupling.
    sst, _, _ = load_year(root, 'sst', y)
    ci, _, _ = load_year(root, 'ci', y)
    if tsrc is not None:
        out['SO clear-sky SW [W/m2]'] = zband(tsrc.mean(axis=0), lat, *SO)
    if sst is not None:
        out['SO SST [C]'] = zband_ocean(sst.mean(axis=0), lat, *SO, lsm) - 273.15
    if ci is not None:
        out['SO sea ice [%]'] = zband_ocean(ci.mean(axis=0), lat, *SO, lsm) * 100.0
    if stl2 is not None:
        out['Siberia DJF soil [C]'] = boxmean(stl2[d].mean(axis=0), lat, lon, SIB,
                                              lsm) - 273.15
        out['Siberia JJA soil [C]'] = boxmean(stl2[j].mean(axis=0), lat, lon, SIB,
                                              lsm) - 273.15
    return out


print(__doc__)
print('=' * 100)

avail = {}
for tag, root in PAIR:
    ys = years_avail(root, '2t')
    avail[tag] = ys
    print(f'  {tag:9s} {len(ys):3d} years on disk: {ys[0]}-{ys[-1]}' if ys
          else f'  {tag:9s} NO OUTPUT')
usable = [y for y in CLEAN if all(y in avail[t] for t, _ in PAIR)]
print(f'\n  matched clean window: {len(usable)} years '
      f'({usable[0]}-{usable[-1]}), excluding {DIRTY[0]}-{DIRTY[-1]}')
if len(usable) < 30:
    print('  WARNING: below the 30-year campaign minimum')

series = {}
for tag, root in PAIR:
    lsm, lat, lon = load_year(root, 'lsm', usable[0])
    if lsm is not None and lsm.ndim == 3:
        lsm = lsm[0]
    rows = []
    for y in usable:
        m = metrics_for_year(root, y, lsm, lat, lon)
        if m:
            rows.append(m)
    series[tag] = rows
    print(f'  {tag:9s} scored {len(rows)} of {len(usable)} years')

keys = list(series[PAIR[0][0]][0].keys())

print('\n' + '=' * 100)
print(f'\nDETECTION THRESHOLDS, from 11E interannual scatter over the same '
      f'{len(usable)} years:  1.96*sd*sqrt(2/n)\n')
base_tag = PAIR[0][0]
thr = {}
for k in keys:
    v = np.array([r[k] for r in series[base_tag]])
    thr[k] = 1.96 * v.std(ddof=1) * np.sqrt(2.0 / len(v))
    print(f'  {k:24s} sd {v.std(ddof=1):8.4f}   threshold +-{thr[k]:.4f}')

print('\n' + '=' * 100)
print(f'\n11G minus 11E, {len(usable)} matched years.  * = resolved\n')
print(f'  {"metric":24s} {"11E":>10s} {"11G":>10s} {"diff":>11s} {"thr":>9s}')
verdict = {}
for k in keys:
    a = np.mean([r[k] for r in series[PAIR[0][0]]])
    b = np.mean([r[k] for r in series[PAIR[1][0]]])
    d = b - a
    sig = abs(d) > thr[k]
    verdict[k] = (a, b, d, sig)
    print(f'  {k:24s} {a:10.3f} {b:10.3f} {d:+10.3f}{"*" if sig else " "} '
          f'{thr[k]:9.3f}')

try:
    with xr.open_dataset(CERESF) as cds:
        clat = cds['lat'].values
        csel = (clat >= SO[0]) & (clat < SO[1])
        cw = np.cos(np.deg2rad(clat[csel]))
        so_obs = float(np.average(
            cds['toa_cre_sw_clim'].values.mean(axis=0)[csel, :].mean(axis=1), weights=cw))
    a, b, _, _ = verdict['SO SW CRE [W/m2]']
    print(f'\n  CERES SO SW CRE {so_obs:.2f}:  11E is {a - so_obs:+.2f} from it, '
          f'11G {b - so_obs:+.2f}')
except Exception as exc:
    print(f'\n  (CERES anchor unavailable: {exc})')

print('\n' + '=' * 100)
print('\nPRE-REGISTERED READING\n')
toa = verdict['net TOA [W/m2]']
so = verdict['SO SW CRE [W/m2]']
jja = verdict['Siberia JJA T2m [C]']
djf = verdict['Siberia DJF T2m [C]']
print(f'  net TOA   {toa[0]:+.3f} -> {toa[1]:+.3f}  ({toa[2]:+.3f}, '
      f'{"resolved" if toa[3] else "NOT resolved"});  AMIP predicted -0.39')
print(f'  SO SW CRE {so[0]:+.3f} -> {so[1]:+.3f}  ({so[2]:+.3f}, '
      f'{"resolved" if so[3] else "NOT resolved"});  AMIP -4.67, '
      f'coupled 1.40x would be about -6.5')
print(f'  Siberia JJA {jja[2]:+.3f} ({"resolved" if jja[3] else "ns"}), '
      f'DJF {djf[2]:+.3f} ({"resolved" if djf[3] else "ns"});  AMIP cost +1.03 / +0.05')
print()
AMIP_SO_CRE, AMIP_SO_NET = -4.704, -3.813   # S4 minus control, evaluate_L.out
net = verdict.get('SO net TOA [W/m2]')
print(f'\n  AMIP anchors (S4 minus control, 48 yr): SO SW CRE {AMIP_SO_CRE:+.3f}, '
      f'SO net TOA {AMIP_SO_NET:+.3f}')
print(f'  coupled / AMIP ratio:  CRE {so[2] / AMIP_SO_CRE:.2f}x, '
      f'net TOA {net[2] / AMIP_SO_NET:.2f}x   (DMS gave 1.40x on net TOA)')
for k in ('SO clear-sky SW [W/m2]', 'SO SST [C]', 'SO sea ice [%]'):
    if k in verdict:
        a, b, d, sg = verdict[k]
        print(f'  {k:24s} {d:+8.3f} ({"resolved" if sg else "ns"})')
print()
if so[3] and so[2] < 0 and toa[2] < 0:
    print('  -> S4 TRANSFERS in sign. Read the ratio and the surface terms above')
    print('     before calling the magnitude: CRE is the part that does NOT survive')
    print('     coupling, so a small CRE ratio is expected, not disqualifying.')
elif not so[3]:
    print('  -> NOT RESOLVED at this length. The lever is not visible above the coupled')
    print('     interannual scatter; do not adopt or reject on this.')
else:
    print('  -> MOVES, BUT NOT AS PREDICTED. Check the sign against AMIP before reading')
    print('     anything into the Siberian columns.')
print()
