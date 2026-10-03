"""Hemispheric sea-ice VOLUME per year, from FESOM's own integrals (sivoln/sivols).

WHY VOLUME.  Extent is geometrically saturated in this configuration and hid an entire
Arctic runaway for 40 model years: the 1950 albedo step showed almost nothing in extent
while volume climbed 24 -> 48 x10^3 km3.  Volume is the state variable that says where
the pack is GOING.  Extent is reported beside it only so the two can be seen to disagree.

sivoln / sivols are FESOM diagnostics, integrated on the native mesh inside the model, so
this needs no mesh file and cannot pick up the wrong one (a real hazard here: CORE3 comes
in a 211567-node and a 220509-node flavour).

Reference values, NH volume x10^3 km3: PIOMAS 1979-2023 gives about 28.7 in April and
13.0 in September.  Those are satellite-era, so a PI run is expected ABOVE them, by an
amount nobody has pinned down; treat them as a floor, not a target.

Usage:  ARMS=PICAL_ccnice,PICAL_momixoff Y0=1941 Y1=1969 \
          python3 scripts/analysis/ice_volume_series.py
"""
import os, glob
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
Y0, Y1 = int(os.environ.get('Y0', 1941)), int(os.environ.get('Y1', 1969))
ARMS = [a.strip() for a in os.environ.get('ARMS', 'PICAL_ccnice').split(',')]

def series(arm, var, y):
    p = f'{R}/{arm}/outdata/fesom/{var}.fesom.{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p, decode_times=False) as d:
        a = np.squeeze(np.asarray(d[var].values, float))
    return a if a.size == 12 else None      # monthly only; skip odd cadences

def trend(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    k = np.isfinite(y)
    if k.sum() < 4: return np.nan, np.nan
    x, y = x[k], y[k]; n = len(x)
    b, a = np.polyfit(x, y, 1)
    r = y - (a + b*x)
    se = np.sqrt((r**2).sum()/(n-2)/((x-x.mean())**2).sum())
    return b*10.0, 1.96*se*10.0

print(__doc__.split('\n')[0]); print(f'window {Y0}-{Y1}, arms {ARMS}')
print('units 10^3 km3 = 10^12 m3\n')
for arm in ARMS:
    yrs, nh, sh = [], [], []
    for y in range(Y0, Y1+1):
        n, s = series(arm, 'sivoln', y), series(arm, 'sivols', y)
        if n is None or s is None: continue
        # sivoln/sivols carry units '1e9 m3', so 10^3 km3 = 10^12 m3 is /1e3, NOT /1e12.
        yrs.append(y); nh.append(n/1e3); sh.append(s/1e3)
    if not yrs: print(f'{arm}: no monthly sivoln/sivols in window\n'); continue
    NH, SH = np.array(nh), np.array(sh)
    rows = [('NH annual', NH.mean(1)), ('NH Apr (max)', NH[:, 3]), ('NH Sep (min)', NH[:, 8]),
            ('SH annual', SH.mean(1)), ('SH Sep (max)', SH[:, 8]), ('SH Feb (min)', SH[:, 1])]
    print(f'=== {arm}  ({yrs[0]}-{yrs[-1]}, {len(yrs)} yr) ===')
    print(f'{"metric":<15}{"first5":>9}{"last5":>9}{"delta":>9}'
          f'{"trend/dec":>12}{"95% CI":>9}   verdict')
    for name, v in rows:
        b, ci = trend(yrs, v)
        sig = 'RUNNING' if np.isfinite(ci) and abs(b) > ci else 'settled'
        print(f'{name:<15}{v[:5].mean():>9.2f}{v[-5:].mean():>9.2f}'
              f'{v[-5:].mean()-v[:5].mean():>9.2f}{b:>12.3f}{ci:>9.3f}   {sig}')
    print(f'  PIOMAS floor: NH Apr 28.7, NH Sep 13.0 (satellite era; PI expected above)\n')
