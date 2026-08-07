"""Offline (SWEMIN, SCALE) sweep against Rutgers, to design P5+ before burning compute.

WHY THIS EXISTS.  Rutgers 24 km says P3 holds +0.259 too much Siberian cover in
September and +0.249 in October, against +0.107/+0.103 for the as-released ramp --
the fitted curve made a known autumn error 2.4x worse.  May is still +0.074 too high.
The two errors have different controls:

  SWEMIN  floors d_c at SWEMIN/rho.  Because it is a floor it binds ONLY at low
          density, i.e. only on fresh autumn snow -- at rho=300 the fitted d_c
          (0.097 m) already exceeds any plausible floor.  It is therefore an
          almost pure AUTUMN control.  It is currently 3.0 kg/m2, chosen as the
          smallest value that beats the as-released numerical safety margin after
          the P1 crash, and never calibrated for its cover role.
  SCALE   multiplies d_c at all densities, so it moves SPRING (where the fitted
          term dominates) more than autumn.  P3 uses 1, P4 uses 3; Rutgers says
          spring wants ~3, soil temperature mildly prefers 1.

So the pair may be separable: raise SWEMIN to fix autumn, raise SCALE to fix spring.
This script tests that offline instead of guessing, which is the step whose absence
produced the P1 crash (full cover on a 0.5 mm pack) and a SWEMIN=5 chosen by hunch.

THE LIMITATION, stated up front because it bounds every number below.  This
reconstructs cover from P3's OWN daily snow state (sd, rsn, cvh) under different
parameters.  It therefore assumes the snowpack does not respond to the cover change
-- but it does: cover sets albedo, which sets melt, which sets SWE.  The feedback is
POSITIVE in autumn (less cover -> more absorbed SW -> less accumulation -> less cover
still), so the true sensitivity is LARGER than shown here and these are LOWER bounds
on the effect of raising SWEMIN.  Use this to rank candidates and pick two, not to
predict the answer.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

RUN = 'amip_P3_scffit'
YEARS = range(1890, 1916)          # 26 yr; enough for a monthly climatology
DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
DOY_MONTH = np.repeat(np.arange(12), DPM)
BASE = dict(dcl=0.014, dch=0.026, md=4.70, bl=1.46, bh=0.40, dcmax=0.30, rhoref=200.0)

# Rutgers 24 km, Siberian box, area-weighted, 1980-2024 (rutgers_vs_model_cover.py)
RUTGERS = np.array([0.999, 1.000, 0.999, 0.956, 0.655, 0.202,
                    0.011, 0.002, 0.067, 0.599, 0.969, 0.997])
# The model box carries a permanent ~0.035 bare fraction (lakes/water in the land
# mask) present in EVERY run including scheme-off, so midwinter cannot reach 1.0.
# Scoring against raw Rutgers would charge every candidate for that; the offset is
# removed so the sweep is judged on the seasons the scheme actually controls.
OFFSET = 0.035
# Score the WHOLE snow season Sep->May, not a hand-picked subset.  A six-month
# score hid the November/December deficit that large SWEMIN causes: raising the
# floor fixes September by delaying saturation, and if you do not score the months
# where saturation SHOULD have happened you never see the price.
SCORE_MONTHS = [8, 9, 10, 11, 0, 1, 2, 3, 4]   # Sep Oct Nov Dec Jan Feb Mar Apr May


def scf(depth, rho, cvh, swemin, scale):
    p = BASE
    r = np.maximum(rho, 50.0) / p['rhoref']
    floor = swemin / np.maximum(rho, 1.0)
    dcl = np.minimum(p['dcmax'], np.maximum(scale * p['dcl'] * r ** p['md'], floor))
    dch = np.minimum(p['dcmax'], np.maximum(scale * p['dch'] * r ** p['md'], floor))
    sl = np.clip((depth / np.maximum(dcl, 1e-9)) ** p['bl'], 0, 1)
    sh = np.clip((depth / np.maximum(dch, 1e-9)) ** p['bh'], 0, 1)
    f = np.clip(cvh, 0.0, 1.0)
    live = (depth > 1e-6) & (depth * rho > 1e-6)
    return np.where(live, np.clip((1 - f) * sl + f * sh, 0, 1), 0.0)


def load(var, y, sel):
    for pat in (f'atm_remapped_1d_{var}_1d_{y}-{y}.nc', f'atm_remapped_1d_{var}_{y}-{y}.nc'):
        f = f'{RT}/{RUN}/outdata/oifs/{pat}'
        if os.path.exists(f):
            with xr.open_dataset(f, decode_times=False) as d:
                a = d[var].values
            if a.shape[0] == 366:
                a = np.delete(a, 59, axis=0)
            return a.reshape(a.shape[0], -1)[:, sel] if a.shape[0] == 365 else None
    return None


with xr.open_dataset(LSMF) as d:
    lsm = d['lsm'].isel(time_counter=0).values
    lat, lon = d['lat'].values, d['lon'].values
LA = np.broadcast_to(lat[:, None], lsm.shape)
LO = np.broadcast_to(lon[None, :], lsm.shape)
sel = np.flatnonzero(((LA >= 55) & (LA <= 75) & (LO >= 60) & (LO <= 180) & (lsm > 0.5)).ravel())
w = np.cos(np.deg2rad(LA)).ravel()[sel]; W = w / w.sum()

print(__doc__.split('THE LIMITATION')[0])
print('=' * 96)

# cache the daily state once; the sweep is then pure arithmetic
DEPTH, RHO, CVH, MSK = [], [], [], []
n = 0
for y in YEARS:
    sd, rsn = load('sd', y, sel), load('rsn', y, sel)
    if sd is None or rsn is None:
        continue
    cvh = load('cvh', y, sel)
    rho = np.maximum(rsn, 1e-6)
    DEPTH.append(sd * 1000.0 / rho); RHO.append(rho)
    CVH.append(np.zeros_like(sd) if cvh is None else cvh)
    n += 1
print(f'  cached {n} yr of daily state from {RUN}\n')
DEPTH = np.concatenate(DEPTH); RHO = np.concatenate(RHO); CVH = np.concatenate(CVH)
MO = np.tile(DOY_MONTH, n)

target = np.clip(RUTGERS - OFFSET, 0, 1)


def monthly(swemin, scale):
    c = scf(DEPTH, RHO, CVH, swemin, scale)
    return np.array([(c[MO == m].mean(axis=0) @ W) for m in range(12)])


print(f'  Rutgers target (minus the {OFFSET:.3f} permanent bare fraction):')
print(f'  {"":11s}' + ''.join(f'{MON[m]:>7s}' for m in SCORE_MONTHS))
print(f'  {"target":11s}' + ''.join(f'{target[m]:7.3f}' for m in SCORE_MONTHS))

print('\n\nSWEEP -- box-mean cover by month, and RMSE over Sep/Oct/Nov/Apr/May/Jun\n')
print(f'  {"SWEMIN":>7s}{"SC":>4s}' + ''.join(f'{MON[m]:>7s}' for m in SCORE_MONTHS)
      + f'{"RMSE":>9s}{"maxerr":>8s}')
rows = []
for swemin in (3.0, 6.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 55.0):
    for scale in (1.0,):
        cov = monthly(swemin, scale)
        err = cov[SCORE_MONTHS] - target[SCORE_MONTHS]
        rms = float(np.sqrt(np.mean(err ** 2)))
        rows.append((rms, swemin, scale, cov))
        print(f'  {swemin:7.1f}{scale:4.0f}' + ''.join(f'{cov[m]:7.3f}' for m in SCORE_MONTHS)
              + f'{rms:9.4f}{err[np.argmax(np.abs(err))]:+8.3f}')

rows.sort()
print(f'\n  BEST BY RMSE: SWEMIN={rows[0][1]:.1f} SCALE={rows[0][2]:.1f}  (RMSE {rows[0][0]:.4f})')
print(f'  P3 AS RUN   : SWEMIN=3.0 SCALE=1.0  (RMSE '
      f'{[r[0] for r in rows if r[1]==3.0 and r[2]==1.0][0]:.4f})')

print("""
  CHOOSING P5/P6.  Pick the best RMSE candidate, and one that trades differently --
  e.g. best-autumn vs best-spring -- rather than two neighbours, so the pair brackets
  the trade-off instead of measuring the same point twice.  Remember the feedback
  caveat: these are LOWER bounds on the autumn sensitivity, so the true optimum
  SWEMIN is probably smaller than the sweep suggests.  And the winter falsifier still
  applies to every candidate: DJF soil must stay at the N1 reference, which SWEMIN
  can only help (a higher floor means less cover on thin snow, never more).""")
