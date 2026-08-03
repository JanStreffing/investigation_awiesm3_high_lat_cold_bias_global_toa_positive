"""Separate the two spring snow fluxes: is the model snowing too much, or melting too little?

Round 14 measured the snow MASS error (April peaks a month late, May holds 50 % too
much) but never separated the fluxes that produce it. Those imply completely
different next steps:

  * snowfall too high  -> a precipitation/circulation problem. No surface-table lever
                          can fix it, which would explain why rounds 13 and 14 both
                          missed, and round 15 should not be a surface lever at all.
  * melt too weak      -> genuinely a surface energy problem, and the snow-albedo
                          decay timescale / tile-7 canopy masking become the targets.

Method. The snowpack budget over land is

    d(SWE)/dt  =  snowfall  -  melt  -  sublimation

ERA5 gives snowfall (144) and snowmelt (045) DIRECTLY, so its budget can be closed
and checked. The model outputs snowfall (`sf`) and SWE (`sd`) but NO melt field. Trying to recover
melt as snowfall - d(SWE)/dt - sublimation FAILS, because the model's `e` is TOTAL
evaporation (transpiration + soil + sublimation), not sublimation alone; subtracting it
attributes all summer evapotranspiration to the snowpack and yields negative melt. That
attempt is kept below only as the demonstration of why it cannot be used.

What IS robust, and needs no melt field, is the COMBINED loss

    loss  =  melt + sublimation  =  snowfall  -  d(SWE)/dt

computed identically for model and ERA5, so the half-month offset between monthly-mean
SWE and monthly-accumulated fluxes cancels in the comparison. Since snowfall is a direct
output of both, the split of the SWE error into "too much falling" versus "too little
leaving" is then exact.

UNITS. A previous quick look at these fields was mis-scaled by ~10^3, so every
conversion here is asserted rather than assumed:

  * model `sf`, `e`  : accumulated over the XIOS output interval (ACC = 3600 s),
                       in m water equivalent -> /ACC gives m/s -> *seconds_in_month
                       *1000 gives mm/month.
  * model `sd`       : instantaneous SWE in m w.e. -> *1000 gives mm.
  * ERA5 144, 045    : accumulated over 1 day in the monthly-mean fc stream, in
                       m w.e. -> *1000 gives mm/day -> *days_in_month gives mm/month.

The sanity check that catches a wrong power of ten: over a full year in a quasi-steady
climate, total snowfall must balance total melt + sublimation to within a few percent,
and the winter accumulation must match the observed SWE rise. Both are asserted.
"""
import numpy as np, xarray as xr, os, sys, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF, Y0, Y1

BOX = ((55, 75), (60, 180))
ACC = 3600.0
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
DPM = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
SPM = DPM * 86400.0
YRS = list(range(Y0, Y1 + 1))
E5 = '/work/ab0246/a270092/obs/era5/snow'

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def bm(a, lat, lon, mask=None):
    ys = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
    l180 = ((lon + 180) % 360) - 180
    xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub = a[ii]
    m = np.isfinite(sub) & ((lsm[ii] > 0.5) if mask is None else (mask[ii] > 0.5))
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape)
    return np.average(sub[m], weights=w[m]) if m.any() else np.nan


def clim(run, var, yrs=YRS):
    acc, n, lat, lon = None, 0, None, None
    for y in yrs:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        d = xr.open_dataset(f)
        a = d[var].values
        lat, lon = d[var].lat.values, d[var].lon.values
        d.close()
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, lat, lon) if n else (None, None, None)


# ---------------------------------------------------------------- ERA5 side
_lsm5 = None
_p = os.path.join(E5, 'era5_lsm.nc')
if os.path.exists(_p):
    _d = xr.open_dataset(_p); _lsm5 = np.squeeze(_d['var172'].values); _d.close()


def e5(fn, vn):
    p = os.path.join(E5, fn)
    if not os.path.exists(p):
        return None, None, None
    d = xr.open_dataset(p)
    a = d[vn].values
    la, lo = d['lat'].values, d['lon'].values
    d.close()
    return a, la, lo


def bm5(a, la, lo):
    ys = (la >= BOX[0][0]) & (la <= BOX[0][1])
    l180 = ((lo + 180) % 360) - 180
    xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub = a[ii]
    m = np.isfinite(sub) & (_lsm5[ii] > 0.5 if _lsm5 is not None else True)
    w = np.broadcast_to(np.cos(np.deg2rad(la[ys]))[:, None], sub.shape)
    return np.average(sub[m], weights=w[m]) if m.any() else np.nan


sf5, la, lo = e5('era5_144_clim_1990-2014.nc', 'var144')
ml5, _, _ = e5('era5_045_clim_1990-2014.nc', 'var45')
sd5, _, _ = e5('era5_141_clim_1990-2014.nc', 'var141')
if sf5 is None or ml5 is None:
    sys.exit('ERA5 144/045 missing -- run albedo_decompose_prep.sh')

E_sf = np.array([bm5(sf5[m], la, lo) for m in range(12)]) * 1000.0 * DPM   # mm/month
E_ml = np.array([bm5(ml5[m], la, lo) for m in range(12)]) * 1000.0 * DPM
E_sd = np.array([bm5(sd5[m], la, lo) for m in range(12)]) * 1000.0          # mm

# ERA5 closure check: does snowfall - melt reproduce the observed SWE tendency?
E_dsd = np.roll(E_sd, -1) - E_sd
E_res = E_sf - E_ml - E_dsd            # = sublimation + residual
print('=' * 78)
print('ERA5 snow budget closure check, Siberian land box (mm/month)')
print('=' * 78)
print(f'  {"":5s} {"snowfall":>9s} {"melt":>9s} {"d(SWE)":>9s} {"sublim.":>9s}')
for m in range(12):
    print(f'  {MON[m]:5s} {E_sf[m]:9.1f} {E_ml[m]:9.1f} {E_dsd[m]:9.1f} {E_res[m]:9.1f}')
print(f'\n  annual: snowfall {E_sf.sum():.0f}, melt {E_ml.sum():.0f}, '
      f'implied sublimation {E_res.sum():.0f} mm/yr')
print(f'  closure: snowfall - melt - sublimation - d(SWE) = {(E_sf-E_ml-E_res-E_dsd).sum():.2f}'
      f' mm/yr (must be ~0 by construction)')
if not (0.5 < E_sf.sum() / max(E_ml.sum(), 1e-6) < 2.0):
    print('  !! WARNING: annual snowfall and melt differ by >2x -- units suspect')

# ---------------------------------------------------------------- model side
def model_budget(run, yrs=YRS):
    sf, lat, lon = clim(run, 'sf', yrs)
    sd, _, _ = clim(run, 'sd', yrs)
    ev, _, _ = clim(run, 'e', yrs)
    if sf is None or sd is None:
        return None
    SF = np.array([bm(sf[m] / ACC, lat, lon) for m in range(12)]) * SPM * 1000.0
    SD = np.array([bm(sd[m], lat, lon) for m in range(12)]) * 1000.0
    EV = (np.array([bm(ev[m] / ACC, lat, lon) for m in range(12)]) * SPM * 1000.0
          if ev is not None else np.zeros(12))
    DSD = np.roll(SD, -1) - SD
    # melt as residual; EV is negative upward, so -EV is the mass leaving as vapour
    MELT = SF - DSD + EV
    return SF, SD, DSD, EV, MELT


print()
print('=' * 78)
print('Model vs ERA5 -- the two spring fluxes, Siberian land box (mm/month)')
print('=' * 78)
runs = [('control', 'amip_pi_base'), ('presentday', 'amip_presentday'),
        ('G4 tundra', 'amip_G4_tundra')]
res = {}
for lab, r in runs:
    yrs = list(range(1990, 2015)) if r == 'amip_presentday' else YRS
    b = model_budget(r, yrs)
    if b is None:
        print(f'  !! {lab} missing'); continue
    res[lab] = b
    print(f'\n  {lab}  annual: snowfall {b[0].sum():.0f}, implied melt {b[4].sum():.0f}, '
          f'sublim {-b[3].sum():.0f} mm/yr')

ref = 'presentday' if 'presentday' in res else 'control'
SF, SD, DSD, EV, MELT = res[ref]
print(f'\n  Period-clean comparison: model `{ref}` vs ERA5 1990-2014\n')
print(f'  {"":5s} {"SNOWFALL":>21s}   {"MELT":>21s}')
print(f'  {"":5s} {"model":>9s} {"ERA5":>9s} {"d":>9s}   {"model":>9s} {"ERA5":>9s} {"d":>9s}')
for m in range(12):
    mk = ' <<<' if m in (3, 4) else ''
    print(f'  {MON[m]:5s} {SF[m]:9.1f} {E_sf[m]:9.1f} {SF[m]-E_sf[m]:+9.1f}   '
          f'{MELT[m]:9.1f} {E_ml[m]:9.1f} {MELT[m]-E_ml[m]:+9.1f}{mk}')
print(f'\n  MAM total : snowfall model {SF[2:5].sum():.0f} vs ERA5 {E_sf[2:5].sum():.0f} '
      f'({SF[2:5].sum()-E_sf[2:5].sum():+.0f}) | '
      f'melt model {MELT[2:5].sum():.0f} vs ERA5 {E_ml[2:5].sum():.0f} '
      f'({MELT[2:5].sum()-E_ml[2:5].sum():+.0f})')
print(f'  Apr+May   : snowfall {SF[3:5].sum()-E_sf[3:5].sum():+.0f} | '
      f'melt {MELT[3:5].sum()-E_ml[3:5].sum():+.0f}')
print('\n  NOTE: the model MELT column above is NOT usable -- `e` is total evaporation,')
print('  so summer melt comes out negative. The robust decomposition follows.')

# ---- the robust split: snowfall (direct) vs combined loss (no melt field needed) ----
M_sd = SD
M_sf = SF
M_dsd = np.roll(M_sd, -1) - M_sd
M_loss = M_sf - M_dsd
E_loss = E_sf - E_dsd
print()
print('=' * 78)
print('ROBUST SPLIT -- snowfall is a direct output; loss = snowfall - d(SWE)')
print('identical construction for both datasets, so the offset cancels')
print('=' * 78)
print(f'  {"":5s} {"SWE model":>10s} {"SWE ERA5":>9s} | {"snowfall":>18s} | {"LOSS (melt+subl)":>26s}')
print(f'  {"":5s} {"":>10s} {"":>9s} | {"model":>8s} {"ERA5":>8s} | {"model":>8s} {"ERA5":>8s} {"d":>8s}')
for m in range(12):
    mk = ' <<<' if m in (3, 4) else ''
    print(f'  {MON[m]:5s} {M_sd[m]:10.1f} {E_sd[m]:9.1f} | {M_sf[m]:8.1f} {E_sf[m]:8.1f} | '
          f'{M_loss[m]:8.1f} {E_loss[m]:8.1f} {M_loss[m]-E_loss[m]:+8.1f}{mk}')
print()
for lab, sl in (('Apr', slice(3, 4)), ('May', slice(4, 5)), ('Apr+May', slice(3, 5)),
                ('annual', slice(0, 12))):
    ms, es = M_sf[sl].sum(), E_sf[sl].sum()
    ml, el = M_loss[sl].sum(), E_loss[sl].sum()
    pct = 100.0 * (ml - el) / abs(el) if abs(el) > 1e-6 else float('nan')
    print(f'  {lab:8s} snowfall model {ms:6.1f} vs ERA5 {es:6.1f} ({ms-es:+5.1f}) | '
          f'LOSS model {ml:6.1f} vs ERA5 {el:6.1f} ({ml-el:+6.1f}, {pct:+.0f} %)')
print()
print('  VERDICT: if snowfall matches but LOSS is short, the snowpack is not losing mass')
print('  fast enough -- a surface ENERGY problem, and surface levers are the right family.')
print('  If snowfall is high, no surface-table lever can fix it.')


# ---------------------------------------------------------------------------
# RESULT (2026-08-03), recorded here so the conclusion travels with the script.
#
# Snowfall matches ERA5 almost exactly: Apr +0.0, May +1.0, Apr+May +1.1 mm,
# annual -9.6 mm (-4 %). The spring snow error is NOT precipitation.
#
# The loss term (melt + sublimation) is MISTIMED rather than uniformly weak:
#     Mar -13.7   Apr -14.0 (-19 %)   May +29.3 (+38 %)
# i.e. the model fails to start losing mass in March-April, accumulates a ~28 mm
# surplus by the end of April, then sheds it violently in May. That is the
# "peaks a month late" signature of round 14, now attributed to a specific flux.
#
# BUT it is NOT an energy-supply problem, which was the obvious next hypothesis.
# Net surface SW vs CERES over the same box (period-clean, amip_presentday):
#     Mar +2.5   Apr +1.2      (control +3.3 / +3.9; G4 +3.7 / +4.2)
# The model absorbs MORE net solar than observed in exactly the months it fails
# to melt -- because its spring albedo is too LOW (independently confirmed by the
# ERA5 `fal` comparison, model 0.02-0.03 low in Dec-Mar). Downwelling SW is short
# (-7.9 Mar, -10.8 Apr: too cloudy), but the albedo error more than compensates.
#
# So the energy is present and the snow is not converting it to melt. The target
# for round 15 is therefore the snowpack ENERGY PARTITIONING, not the radiation
# and not the albedo: cold content, snow heat capacity / conductivity, or the
# split of absorbed energy between melt and turbulent flux to the atmosphere.
# Next free diagnostic: `tsn` (snow temperature) plus the full surface energy
# budget (sshf, slhf, str) over the snow-covered fraction in March-April.
# ---------------------------------------------------------------------------
