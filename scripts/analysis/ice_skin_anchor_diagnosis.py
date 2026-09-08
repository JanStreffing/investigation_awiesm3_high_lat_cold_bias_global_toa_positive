#!/usr/bin/env python3
"""Test the ist_ref anchor hypothesis against 11X output (core3 mesh).

gen_forcing_couple.F90 stores send field 4, which is ice_temp*a_ice, as
ice%atmcoupl%ist_ref.  ice_thermo_cpl.F90 then uses it as `tref`, the
linearisation anchor of the implicit skin solve.  The anchor therefore carries
a systematic (1-a_ice)*ist cold offset every coupling step.  That offset enters
the solve with gain zlam/(zcpdte + con/zsniced + zlam), which tends to 1 where
snow is thick, so the prediction is that diverged points are thick-snow points.
"""
import glob, numpy as np, netCDF4 as nc

FILL = 1e30
RUN = '/work/bb1469/a270092/runtime/awiesm3-v3.4/11X/run_*/work'
N_CORE3 = 220509


def load(var, year, stream=None):
    """stream lets a variable be read from a differently-named file, e.g. m_ice
    is monthly in m_ice.fesom but daily in m_ice_day.fesom."""
    f = sorted(glob.glob(f'{RUN}/{stream or var}.fesom_{year}-{year}.nc'))
    a = np.asarray(nc.Dataset(f[0]).variables[var][:], dtype=np.float64)
    assert a.shape[1] == N_CORE3, f'{var}: {a.shape[1]} nodes, expected core3 {N_CORE3}'
    return np.where(np.abs(a) > FILL, np.nan, a)


for year in ('1351', '1350'):
    ist, aice, snow = (load(v, year) for v in ('ist', 'a_ice', 'm_snow'))
    ice = aice > 0.01
    cold = ice & (ist < 200.0)
    warm = ice & (ist > 240.0)
    print(f'=== {year} ===  ice point-days {ice.sum():,}   cold(<200K) {cold.sum():,}')
    for name, m in (('cold <200K', cold), ('normal >240K', warm)):
        if not m.any():
            continue
        sn = snow[m] / np.maximum(aice[m], 1e-8)   # grid-mean -> per-ice depth
        print(f'  {name:13s} a_ice med={np.nanmedian(aice[m]):.3f}'
              f'  snow/ice med={np.nanmedian(sn):.3f} m'
              f'  ist med={np.nanmedian(ist[m]):.1f} K')

    # The anchor offset (1-a_ice)*ist should predict the next day's cooling.
    cur, nxt = ist[:-1], ist[1:]
    offset = (1.0 - aice[:-1]) * ist[:-1]
    m = (aice[:-1] > 0.5) & (cur > 200) & np.isfinite(cur) & np.isfinite(nxt)
    drop = (cur - nxt)[m]
    off = offset[m]
    print(f'  anchor offset (1-a)*ist: med={np.median(off):.2f} K'
          f'  p90={np.percentile(off, 90):.2f} K'
          f'  corr with next-day cooling r={np.corrcoef(off, drop)[0, 1]:+.3f}'
          f'  (n={m.sum():,})')

# --- Snow insulation profile -------------------------------------------------
# Thick snow shrinks con/zsniced, so the ice interior loses its hold on the skin
# and whatever error sits in a2ihf sets the temperature.  If that is the enabling
# condition, ist should fall monotonically with snow depth.
print('\n=== ist vs per-ice snow depth (1350, a_ice>0.5) ===')
ist, aice, snow = (load(v, '1350') for v in ('ist', 'a_ice', 'm_snow'))
m = (aice > 0.5) & np.isfinite(ist) & np.isfinite(snow)
sn = snow[m] / aice[m]
t = ist[m]
edges = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 99]
for lo, hi in zip(edges[:-1], edges[1:]):
    b = (sn >= lo) & (sn < hi)
    if b.sum() < 100:
        continue
    print(f'  snow {lo:5.2f}-{hi:5.2f} m  n={b.sum():10,}  ist med={np.median(t[b]):6.1f} K'
          f'  p1={np.percentile(t[b],1):6.1f}  frac<200K={100*(t[b]<200).mean():5.2f}%')

# --- The right statistic: spurious conductance vs real conduction -------------
# With tref = a_ice*t (the weighting bug), the steady state of ice_surftemp is
#   t* = [(con/zsniced)*TFrezs + a2ihf] / [con/zsniced + zlam*(1-a_ice)]
# so the weighting acts as an extra conductance zlam*(1-a_ice) toward 0 K.  The
# controlling parameter is the RATIO of that to the real conduction, not the
# offset (1-a)*ist.  Prediction: diverged points have R >~ 1, healthy points R << 1.
CON, CONSN, EMISS, SIGMA, TFREZS = 2.1656, 0.31, 0.97, 5.67e-8, 271.35
print('\n=== spurious/real conductance ratio R (1350, a_ice>0.5) ===')
ist, aice, snow = (load(v, '1350') for v in ('ist', 'a_ice', 'm_snow'))
mice = load('m_ice', '1350', stream='m_ice_day')
m = (aice > 0.5) & np.isfinite(ist) & np.isfinite(snow) & np.isfinite(mice)
a, t = aice[m], ist[m]
h, hsn = mice[m] / a, snow[m] / a                    # grid-mean -> per-ice
zsniced = np.maximum(h + (CON / CONSN) * hsn, 1e-3)
kreal = CON / zsniced                                 # real conduction  [W/m2/K]
tref = a * t                                          # the anchor as stored
zlam = 4.0 * EMISS * SIGMA * tref**3 + 16.0
kspur = zlam * (1.0 - a)                              # spurious conductance
R = kspur / kreal
tstar = (kreal * TFREZS) / (kreal + kspur)            # fixed point at a2ihf = 0

for name, sel in (('cold  <200K', t < 200), ('normal >240K', t > 240)):
    print(f'  {name}  R med={np.median(R[sel]):7.3f}  '
          f'k_real={np.median(kreal[sel]):6.3f}  k_spur={np.median(kspur[sel]):6.3f}  '
          f't* med={np.median(tstar[sel]):6.1f} K  ist med={np.median(t[sel]):6.1f} K')
print(f'  corr(log R, ist) = {np.corrcoef(np.log10(np.maximum(R,1e-6)), t)[0,1]:+.3f}')
print(f'  frac<200K by R:  R<0.1 {100*(t[R<0.1]<200).mean():5.2f}%   '
      f'0.1-1 {100*(t[(R>=0.1)&(R<1)]<200).mean():5.2f}%   '
      f'R>1 {100*(t[R>=1]<200).mean():5.2f}%')
