"""Do AMIP runs drift?  And is year 1870 -- the year every screen used -- representative?

TWO WORRIES, BOTH ANSWERABLE FROM DATA ALREADY ON DISK.

1. SPIN-UP.  runs.py evaluates 1872-1915 and discards 1870-71 "for deep-soil spin-up".
   Every one-year screen this campaign has run -- T1/T2/T3, the U bracket, V, the W
   series, the X ladder, Y -- used YEAR 1870, i.e. the first discarded year.  The defence
   is that a screen is a PAIR at the same year from the same initial state, so a common
   transient cancels in the difference.  That is an argument, not a measurement.  If the
   model is still adjusting, the RESPONSE to a lever can differ from its equilibrated
   response even though the transient itself cancels.

2. DRIFT.  These runs do NOT use fixed 1850 forcing.  NCMIPFIXYR sits in the dead
   NAMECECMIP6 block, so the runtime value is -1 and the model reads yearly CMIP GHG:
   1869 in the first leg, 1911 in the last.  CO2 rises roughly 288 -> 301 ppm across the
   record, about +0.23 W/m2 of forcing.  In AMIP the SST cannot absorb that, so it should
   appear directly at the top of atmosphere as a trend -- and a 44-year MEAN of a
   trending series is not a clean pre-industrial estimate.

WHAT THIS PRINTS, for the control and for one lever arm:
  * the yearly series of each metric, so a transient is visible rather than inferred
  * first-5-year mean vs last-5-year mean, and a least-squares trend per decade
  * where year 1870 sits relative to the evaluated 1872-1915 mean, in units of the
    interannual sd -- the direct test of whether the screening year is representative
  * the same for a LEVER's RESPONSE: does the 1870 delta match the 1872-1915 delta?
    That is the question the screens actually depend on, and it is checked against A1a
    and D2a, the two arms the report already validated at one year.

Metrics are radiative because that is what the one-year screens claim to resolve.  The
Siberian T2m series is included as a contrast: it is the slowest-adjusting field and the
one no screen may touch.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import sys
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runs import RT, Y0, Y1

ACC = 3600.0
CTRL = 'amip_pi_base'
ARMS = ['amip_A1_overlap01', 'amip_D2a_inpsea02']   # validated at 1 yr in the report
SCREEN_YEAR = 1870
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/'
        'Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/outdata/oifs/'
        'atm_remapped_1m_lsm_1350-1350.nc')

print(__doc__)
print('=' * 100)


def series(run, y0, y1):
    """Per-year values of each metric; None entries for missing years."""
    out = {}
    lsm = None
    if os.path.exists(LSMF):
        with xr.open_dataset(LSMF, decode_times=False) as d:
            lsm = np.squeeze(d['lsm'].values)
    years = []
    rows = {k: [] for k in ('global net TOA', 'SO SW CRE', 'trop LW CRE', 'Siberia JJA T2m')}
    for y in range(y0, y1 + 1):
        D = f'{RT}/{run}/outdata/oifs'
        need = {v: f'{D}/atm_remapped_1m_{v}_1m_{y}-{y}.nc'
                for v in ('tsr', 'tsrc', 'ttr', 'ttrc', '2t')}
        if not all(os.path.exists(p) for p in need.values()):
            continue
        d = {}
        for v, p in need.items():
            with xr.open_dataset(p, decode_times=False) as ds:
                d[v] = ds[v].values / (ACC if v in ('tsr', 'tsrc', 'ttr', 'ttrc') else 1.0)
                lat, lon = ds['lat'].values, ds['lon'].values
        sw, lw, net = d['tsr'] - d['tsrc'], d['ttr'] - d['ttrc'], d['tsr'] + d['ttr']

        def b(a, lo, hi):
            s = (lat >= lo) & (lat < hi)
            w = np.cos(np.deg2rad(lat[s]))
            return float(np.average(a.mean(axis=0)[s, :].mean(axis=1), weights=w))
        years.append(y)
        rows['global net TOA'].append(b(net, -90, 90))
        rows['SO SW CRE'].append(b(sw, -65, -45))
        rows['trop LW CRE'].append(b(lw, -30, 30))
        # Siberia 55-75N, 60-180E, land, JJA
        t = d['2t'][[5, 6, 7]].mean(axis=0) - 273.15
        sy = (lat >= 55) & (lat < 75)
        sx = (lon >= 60) & (lon <= 180)
        m = np.ones_like(t, dtype=bool)
        m[:] = False
        m[np.ix_(sy, sx)] = True
        if lsm is not None and lsm.shape == t.shape:
            m &= lsm > 0.5
        w2 = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], t.shape)
        rows['Siberia JJA T2m'].append(float(np.average(t[m], weights=w2[m])))
    out['years'] = np.array(years)
    for k in rows:
        out[k] = np.array(rows[k])
    return out


c = series(CTRL, 1870, 1925)
if len(c['years']) == 0:
    raise SystemExit('no control years found')
yrs = c['years']
print(f'control {CTRL}: {len(yrs)} years, {yrs[0]}-{yrs[-1]}; '
      f'evaluation window {Y0}-{Y1}\n')

KEYS = ['global net TOA', 'SO SW CRE', 'trop LW CRE', 'Siberia JJA T2m']
print('1. IS THERE A TREND?  (least squares over the evaluation window)')
print('-' * 100)
print(f'  {"metric":18s} {"first 5 yr":>11s} {"last 5 yr":>11s} {"trend/decade":>13s} '
      f'{"sd(interann)":>13s} {"trend/sd":>9s}')
win = (yrs >= Y0) & (yrs <= Y1)
for k in KEYS:
    v = c[k][win]
    t = yrs[win]
    sl = np.polyfit(t, v, 1)[0] * 10.0
    sd = v.std(ddof=1)
    print(f'  {k:18s} {v[:5].mean():11.3f} {v[-5:].mean():11.3f} {sl:13.3f} '
          f'{sd:13.3f} {sl/sd if sd else 0:9.2f}')

print('\n2. IS THE SPIN-UP VISIBLE?  first years vs the evaluated mean, in sd units')
print('-' * 100)
print(f'  {"metric":18s} ' + ' '.join(f'{y:>9d}' for y in yrs[:5])
      + f' | {"mean " + str(Y0) + "-" + str(Y1):>14s}')
for k in KEYS:
    v = c[k]
    mu, sd = v[win].mean(), v[win].std(ddof=1)
    z = [(v[i] - mu) / sd for i in range(min(5, len(v)))]
    print(f'  {k:18s} ' + ' '.join(f'{x:+9.2f}' for x in z) + f' | {mu:14.3f}')
print('  (values are z-scores of each year against the evaluated mean; |z|>2 is unusual)')

print(f'\n3. IS YEAR {SCREEN_YEAR} -- THE SCREENING YEAR -- REPRESENTATIVE?')
print('-' * 100)
i0 = int(np.where(yrs == SCREEN_YEAR)[0][0]) if SCREEN_YEAR in yrs else None
if i0 is None:
    print(f'  year {SCREEN_YEAR} not on disk')
else:
    for k in KEYS:
        v = c[k]
        mu, sd = v[win].mean(), v[win].std(ddof=1)
        print(f'  {k:18s} {SCREEN_YEAR}: {v[i0]:8.3f}   mean {mu:8.3f}   '
              f'z = {(v[i0]-mu)/sd:+.2f}')

print('\n4. DOES A LEVER RESPONSE AT YEAR 1870 MATCH ITS 44-YEAR RESPONSE?')
print('-' * 100)
print(f'  {"arm":22s} {"metric":18s} {"1870 delta":>11s} {"44-yr delta":>12s} {"ratio":>8s}')
for arm in ARMS:
    a = series(arm, 1870, 1925)
    if len(a['years']) == 0:
        print(f'  {arm}: not on disk')
        continue
    aw = (a['years'] >= Y0) & (a['years'] <= Y1)
    ai = int(np.where(a['years'] == SCREEN_YEAR)[0][0]) if SCREEN_YEAR in a['years'] else None
    for k in KEYS:
        if ai is None:
            continue
        d1 = a[k][ai] - c[k][i0]
        d44 = a[k][aw].mean() - c[k][win].mean()
        r = d1 / d44 if abs(d44) > 1e-9 else float('nan')
        print(f'  {arm.replace("amip_",""):22s} {k:18s} {d1:11.3f} {d44:12.3f} {r:8.2f}')
print('\n  A ratio near 1 means the one-year screen reproduces the equilibrated response.')
print('  The report already found 1.02 for A1a and 1.09 for D2a on SO SW CRE; anything')
print('  far from 1, or of the wrong sign, marks a metric the screens cannot be trusted on.')
