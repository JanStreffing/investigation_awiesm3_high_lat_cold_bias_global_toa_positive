#!/usr/bin/env python3
"""Decide H1 for the coupled-slab ice skin divergence.

ice_surftemp solves, at steady state,

    t = [a2ihf + zlam*tref + (con/zsniced)*TFrezs] / (con/zsniced + zlam)

with zlam = 4*eps*sigma*tref**3 + zlam_turb and zlam_turb = 16.0 hardcoded.  The
denominator is 19.9-20.5 W/m2/K over the whole observed range of snow and ice
thickness, so the sensitivity to the atmospheric flux is ~0.050 K per W/m2 and
essentially nothing else in the equation is soft enough to move the skin far.

Reaching the observed ~190 K from a ~250 K anchor therefore requires

    a2ihf ~ 190*19.87 - 19.44*250 - 0.433*271.35 ~ -1200 W/m2

H1 (a2ihf is kilowatt-scale, a Bug IV recurrence) predicts exactly that.  If
instead a2ihf at the diverged points is of ordinary magnitude, H1 is dead and
the fault is inside the solve rather than in what the atmosphere delivers.

Two earlier candidate mechanisms were tested against 11X and both failed with
the wrong sign, so they are not re-tested here:
  - the (1-a_ice) offset in ist_ref correlates with next-day cooling at r=-0.003
  - the spurious/real conductance ratio is LARGER at healthy points (1.69 vs 0.76)
"""
import glob, sys, numpy as np, netCDF4 as nc

FILL = 1e30
CON, CONSN, EMISS, SIGMA, TFREZS, ZLAM_TURB = 2.1656, 0.31, 0.97, 5.67e-8, 271.35, 16.0
RUN = '/work/bb1469/a270092/runtime/awiesm3-v3.4/11Xdbg/run_*/work'
N_CORE3 = 220509


def load(stream, var):
    f = sorted(glob.glob(f'{RUN}/{stream}.fesom_*.nc'))
    if not f:
        sys.exit(f'no files for {stream}; has the run written anything yet?')
    out = []
    for p in f:
        d = nc.Dataset(p)
        if var not in d.variables or d.variables[var].shape[0] == 0:
            continue
        out.append(np.asarray(d.variables[var][:], dtype=np.float64))
    if not out:
        sys.exit(f'{stream}: files exist but hold no records (sync_freq not working?)')
    a = np.concatenate(out, axis=0)
    assert a.shape[1] == N_CORE3, f'{stream}: {a.shape[1]} nodes, expected core3 {N_CORE3}'
    return np.where(np.abs(a) > FILL, np.nan, a)


a2ihf = load('a2ihf', 'a2ihf')
istref = load('istref', 'istref')
ist = load('dbg_ist', 'ist')
aice = load('dbg_a_ice', 'a_ice')
snow = load('dbg_m_snow', 'm_snow')
n = min(x.shape[0] for x in (a2ihf, istref, ist, aice, snow))
a2ihf, istref, ist, aice, snow = (x[:n] for x in (a2ihf, istref, ist, aice, snow))
print(f'{n} days of synced daily output\n')

ok = (aice > 0.5) & np.isfinite(a2ihf) & np.isfinite(ist)
cold, warm = ok & (ist < 200), ok & (ist > 240)
print(f'{"population":14s} {"n":>10s} {"a2ihf med":>11s} {"a2ihf p1":>10s} '
      f'{"a2ihf min":>11s} {"ist med":>8s} {"istref med":>11s}')
for name, m in (('cold  <200K', cold), ('normal >240K', warm)):
    if not m.any():
        print(f'  {name:12s} (none)')
        continue
    print(f'  {name:12s} {m.sum():10,} {np.nanmedian(a2ihf[m]):11.1f} '
          f'{np.nanpercentile(a2ihf[m],1):10.1f} {np.nanmin(a2ihf[m]):11.1f} '
          f'{np.nanmedian(ist[m]):8.1f} {np.nanmedian(istref[m]):11.1f}')

# Does the solve's own steady state reproduce the observed ist from the observed
# a2ihf?  If it does, the flux is the whole story and the solve is innocent.
h_sn = snow / np.maximum(aice, 1e-8)
zsniced = np.maximum(0.5 + (CON / CONSN) * h_sn, 1e-3)   # h unavailable here; 0.5 m nominal
kreal = CON / zsniced
zlam = 4 * EMISS * SIGMA * np.maximum(istref, 1.0)**3 + ZLAM_TURB
pred = (a2ihf + zlam * istref + kreal * TFREZS) / (kreal + zlam)
for name, m in (('cold  <200K', cold), ('normal >240K', warm)):
    if m.any():
        print(f'  {name:12s} predicted ist from a2ihf: {np.nanmedian(pred[m]):7.1f} K '
              f'(observed {np.nanmedian(ist[m]):.1f} K)')

print(f'\nistref vs ist*a_ice (weighting check): '
      f'median ratio istref/ist = {np.nanmedian((istref/np.maximum(ist,1))[ok]):.4f}, '
      f'median a_ice = {np.nanmedian(aice[ok]):.4f}')
med = np.nanmedian(a2ihf[cold]) if cold.any() else np.nan
print(f'\nVERDICT: a2ihf at diverged points = {med:.1f} W/m2 -> '
      f'{"H1 CONFIRMED (kilowatt-scale)" if med < -500 else "H1 DEAD, fault is inside the solve"}')
