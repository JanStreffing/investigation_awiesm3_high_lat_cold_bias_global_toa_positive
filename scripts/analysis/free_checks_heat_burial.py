"""Four checks that could redirect the whole "reduce heat burial" plan before any lever.

  1  WHERE DOES THE HEAT ENTER?   Surface heat flux per band against the heat that band
     actually gained.  The residual is lateral convergence.  We measured where the heat
     LANDS; we never measured where it crosses the surface.
     NOTE ON SIGN: FESOM's `fh` is positive UPWARD (out of the ocean).  Its global mean
     over 1380-89 is -1.535 W/m2 against a storage rate of +1.455 from thetaoga, so the
     budget closes to 0.08 W/m2 only with that convention.  An earlier pass took it as
     positive-downward and produced "lateral" terms with the wrong sign throughout.
  2  IS THE DRIFT DECELERATING?   thetaoga over all 40 years, decade by decade.  Steady
     means a bias; decaying means we would be tuning against spin-up.
  3  IS SALT BEING CONSERVED?     which_ale='linfs' applies freshwater as a VIRTUAL SALT
     FLUX.  If the global salt inventory drifts, the surface salinity trends at the
     outcrops (+0.30 psu at 45-60N over 30 yr) are a numerical artefact -- and they, not
     the mixing coefficients, would be what deepens the winter mixed layer there.
  4  DID THE WINDS CHANGE?        The tropical thermocline sank 42 m.  That is Ekman
     pumping.  If the stress drifted within the run, the heave is dynamic.

  5  BONUS, and the sharpest of them: opottempdiff vs opottemprmadvect, the CMOR
     decomposition of the column heat tendency into parameterised dianeutral mixing and
     advection.  This measures directly what the isopycnal analysis could only infer.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

D = '/work/bb1469/a270092/runtime/awiesm3-v3.4/11X/outdata/fesom'
MESHF = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
YEARS = list(range(1350, 1390))
RHO_CP = 1026.0 * 3996.0
SEC_Y = 365.0 * 86400.0
BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45), ('45-30S', -45, -30),
         ('30S-30N', -30, 30), ('30-45N', 30, 45), ('45-60N', 45, 60), ('60-90N', 60, 90)]

m = xr.open_dataset(MESHF, decode_times=False)
lat = m['lat'].values.astype('f8')
area = m['cell_area'].values.astype('f8')
tri = m['triag_nodes'].values.astype(int) - 1
elat = lat[tri].mean(axis=1)
assert lat.size == 220509, 'MESH GUARD'


def series(var):
    out = []
    for y in YEARS:
        f = f'{D}/{var}.fesom.{y}.nc'
        if not os.path.exists(f):
            out.append(np.nan); continue
        with xr.open_dataset(f, decode_times=False) as d:
            out.append(float(d[var].values.mean()))
    return np.array(out)


def field_mean(var, years):
    acc, n = None, 0
    for y in years:
        f = f'{D}/{var}.fesom.{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values.astype('f8').mean(axis=0)
        acc = a if acc is None else acc + a
        n += 1
    return acc / max(n, 1), n


def bandmean(f, la, ar):
    out = []
    for nm, lo, hi in BANDS:
        s = (la >= lo) & (la < hi) & np.isfinite(f)
        out.append(np.sum(f[s] * ar[s]) / np.sum(ar[s]) if s.any() else np.nan)
    return np.array(out)


# =============================================================== 2  deceleration
print('=' * 100)
print('CHECK 2  IS THE DRIFT DECELERATING?   global mean potential temperature, 11X')
print('=' * 100)
th = series('thetaoga')
vol = series('volo')
ga = area.sum()
print(f'  ocean volume {np.nanmean(vol):.4e} m3, surface area used {ga:.4e} m2')
print(f'{"decade":<14}{"thetaoga [degC]":>18}{"d/dt [K/decade]":>18}{"implied [W m-2]":>18}')
prev = None
for i in range(0, 40, 10):
    dec = th[i:i + 10]
    yy = YEARS[i]
    mean = np.nanmean(dec)
    if prev is not None:
        dK = mean - prev
        w = dK * RHO_CP * np.nanmean(vol) / ga / (10 * SEC_Y)
        print(f'{yy}-{yy+9:<9}{mean:>18.5f}{dK:>+18.5f}{w:>+18.3f}')
    else:
        print(f'{yy}-{yy+9:<9}{mean:>18.5f}{"--":>18}{"--":>18}')
    prev = mean
tot = (np.nanmean(th[30:40]) - np.nanmean(th[0:10])) * RHO_CP * np.nanmean(vol) / ga / (30 * SEC_Y)
print(f'  1350-59 -> 1380-89 overall: {tot:+.3f} W m-2 of ocean area')

# =============================================================== 3  salt budget
print('\n' + '=' * 100)
print('CHECK 3  IS SALT CONSERVED?   linfs applies freshwater as a virtual salt flux')
print('=' * 100)
so = series('soga')
print(f'{"decade":<14}{"soga [psu]":>14}{"drift [psu/decade]":>22}')
prev = None
for i in range(0, 40, 10):
    mean = np.nanmean(so[i:i + 10])
    d = f'{mean - prev:+.6f}' if prev is not None else '--'
    print(f'{YEARS[i]}-{YEARS[i]+9:<9}{mean:>14.5f}{d:>22}')
    prev = mean
print(f'  total soga drift 1350-59 -> 1380-89: {np.nanmean(so[30:40])-np.nanmean(so[0:10]):+.6f} psu')
vs, _ = field_mean('virtsalt', YEARS[-10:])
rs, _ = field_mean('relaxsalt', YEARS[-10:])
print(f'  global mean virtsalt {np.nansum(vs*area)/ga:+.4e} m/s psu   '
      f'implied soga tendency {np.nansum(vs*area)*SEC_Y/np.nanmean(vol):+.6f} psu/yr')
print(f'  global mean relaxsalt {np.nansum(rs*area)/ga:+.4e} m/s psu  '
      f'(surf_relax_s is 0, so this should be 0)')
print(f'\n{"band":<12}{"virtsalt [m/s psu]":>22}{"as psu/decade in top 100 m":>30}')
bvs = bandmean(vs, lat, area)
for i, (nm, lo, hi) in enumerate(BANDS):
    print(f'{nm:<12}{bvs[i]:>22.4e}{bvs[i]*10*SEC_Y/100.0:>30.4f}')

# =============================================================== 1  surface vs gained
print('\n' + '=' * 100)
print('CHECK 1  WHERE DOES THE HEAT ENTER?   surface flux vs heat actually gained')
print('=' * 100)
fh_e, _ = field_mean('fh', YEARS[:10])
fh_l, _ = field_mean('fh', YEARS[-10:])
d = np.load('data/ohc_drift_11X.npz')
dT, dzl = d['dT'], d['dz']
wet = np.isfinite(dT)
volc = area[:, None] * dzl[None, :] * wet
print('  fh is positive UPWARD, so "into ocean" = -fh')
print(f'{"band":<12}{"into ocean":>12}{"lateral":>11}{"kept":>10}')
rows = []
for nm, lo, hi in BANDS:
    s = (lat >= lo) & (lat < hi)
    ba = area[s & np.any(wet, axis=1)].sum()
    q = np.nansum(dT[s] * volc[s]) * RHO_CP / ba / (30 * SEC_Y)
    into = -np.nansum(fh_l[s] * area[s]) / area[s].sum()
    print(f'{nm:<12}{into:>+12.2f}{q - into:>+11.2f}{q:>+10.3f}')
    rows.append((nm, into, q - into, q))
print(f'{"GLOBAL":<12}{-np.nansum(fh_l*area)/ga:>+12.2f}')
print('  "lateral" = kept minus taken in through the surface, i.e. ocean transport')

# =============================================================== 5  diff vs advect
print('\n' + '=' * 100)
print('CHECK 5  WHAT MOVES THE HEAT?  column-integrated tendency, late decade [W m-2]')
print('=' * 100)
di, _ = field_mean('opottempdiff', YEARS[-10:])
ad, _ = field_mean('opottemprmadvect', YEARS[-10:])
te, _ = field_mean('opottemptend', YEARS[-10:])
print(f'{"band":<12}{"diffusion":>13}{"advection":>13}{"net tendency":>15}')
for i, (nm, lo, hi) in enumerate(BANDS):
    s = (lat >= lo) & (lat < hi)
    w = area[s]
    print(f'{nm:<12}{np.nansum(di[s]*w)/w.sum():>+13.3f}'
          f'{np.nansum(ad[s]*w)/w.sum():>+13.3f}{np.nansum(te[s]*w)/w.sum():>+15.3f}')
print(f'{"GLOBAL":<12}{np.nansum(di*area)/ga:>+13.3f}{np.nansum(ad*area)/ga:>+13.3f}'
      f'{np.nansum(te*area)/ga:>+15.3f}')

np.savez('data/free_checks_11X.npz', theta=th, soga=so, volo=vol,
         budget=np.array([[r[1], r[2], r[3]] for r in rows]),
         bands=np.array([b[0] for b in BANDS]), years=np.array(YEARS))

# =============================================================== 4  winds
print('\n' + '=' * 100)
print('CHECK 4  DID THE WINDS DRIFT?   zonal surface stress to ocean [N m-2]')
print('=' * 100)
tx_e, _ = field_mean('tx_sur', YEARS[:10])
tx_l, _ = field_mean('tx_sur', YEARS[-10:])
ea = np.ones_like(elat)
print(f'{"band":<12}{"early":>12}{"late":>12}{"change":>12}{"% change":>11}')
be = bandmean(tx_e, elat, ea)
bl = bandmean(tx_l, elat, ea)
for i, (nm, lo, hi) in enumerate(BANDS):
    pc = 100 * (bl[i] - be[i]) / abs(be[i]) if abs(be[i]) > 1e-6 else np.nan
    print(f'{nm:<12}{be[i]:>+12.4f}{bl[i]:>+12.4f}{bl[i]-be[i]:>+12.5f}{pc:>+11.2f}')
