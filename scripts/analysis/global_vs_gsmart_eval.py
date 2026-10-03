#!/usr/bin/env python3
"""Evaluate OASIS CONSERV GLOBAL against GSMART for A_Q_ice.

A_Q_ice carries CONSERV/GSMART, which resolves to GLBPOS: a MULTIPLICATIVE
correction zlagr = av_sums/av_sumd applied to every destination point
(mod_oasis_advance.F90). A_Q_ice is sign-changing - positive over summer-hemisphere
ice, negative over winter - so its global integral passes through zero twice a
year and takes the denominator with it. GLOBAL instead subtracts the area-mean
residual, which is bounded by the remap's conservation error and cannot blow up.

Both runs are the same code (skin taken from OpenIFS, recovery fix, Waitgroup fix)
and the same year from the same restart. They differ in the CONSERV method on
A_Q_ice. 11Xg also drops a stray CONSERV from the temperature anchor, but ist
moved only 0.13 K between the runs, so the mass-budget difference is attributable
to A_Q_ice.

Prediction of the mechanism: a multiplicative blow-up is tied to the global
integral crossing zero, so the difference should CONCENTRATE at the seasonal
transitions rather than spread evenly through the year.
"""
import glob, sys, numpy as np, netCDF4 as nc

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
RUNS = {'GSMART (E4)': f'{R}/11Xdbg/run_*/work', 'GLOBAL (11Xg)': f'{R}/11Xg/run_*/work'}


def series(base, stream):
    """Return the time series of a scalar CMOR-style diagnostic.

    The variable carries the stream's own name; select it explicitly rather than
    by elimination, or 'time_bounds' (12, 2) gets picked up instead.
    """
    out = []
    for p in sorted(glob.glob(f'{base}/{stream}.fesom_*.nc')):
        d = nc.Dataset(p)
        if stream in d.variables and d.variables[stream].shape[0]:
            out.append(np.asarray(d.variables[stream][:], dtype=np.float64).squeeze())
    if not out:
        return None
    a = np.concatenate([np.atleast_1d(x) for x in out])
    return np.where(np.abs(a) > 1e30, np.nan, a)


def field(base, stream, var):
    """Return a (time, node) field."""
    out = []
    for p in sorted(glob.glob(f'{base}/{stream}.fesom_*.nc')):
        d = nc.Dataset(p)
        if var in d.variables and d.variables[var].shape[0]:
            out.append(np.asarray(d.variables[var][:], dtype=np.float64))
    if not out:
        return None
    a = np.concatenate(out)
    return np.where(np.abs(a) > 1e30, np.nan, a)


print(f'  {"diagnostic":12s} {"run":14s} {"start":>12s} {"end":>12s} {"drift":>10s} {"mean":>12s}')
store = {}
for stream in ('sivoln', 'sivols', 'siarean', 'siareas', 'thetaoga'):
    for name, pat in RUNS.items():
        base = glob.glob(pat)
        if not base:
            continue
        s = series(base[0], stream)
        if s is None or s.size < 2:
            continue
        store[(stream, name)] = s
        d = s[-1] - s[0]
        rel = 100 * d / abs(s[0]) if np.isfinite(s[0]) and s[0] != 0 else np.nan
        print(f'  {stream:12s} {name:14s} {s[0]:12.5g} {s[-1]:12.5g} '
              f'{rel:9.2f}% {np.nanmean(s):12.5g}')
    print()

# Where in the year do the two runs diverge?
print('  Monthly divergence, GLOBAL minus GSMART (relative to the GSMART mean):')
names = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
for stream in ('sivoln', 'sivols'):
    a = store.get((stream, 'GSMART (E4)'))
    b = store.get((stream, 'GLOBAL (11Xg)'))
    if a is None or b is None:
        continue
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    ref = np.nanmean(a)
    row = 100 * (b - a) / ref
    print(f'   {stream}: ' + ' '.join(f'{nm}{v:+6.2f}' for nm, v in zip(names, row)))

# Thermodynamic ice growth is driven by Qatmice = -a2ihf, so it sees the flux
# difference without needing the flux itself on disk.
print()
print('  thdgrice (thermodynamic ice growth), monthly global mean [m/s]:')
ga = field(glob.glob(RUNS['GSMART (E4)'])[0], 'thdgrice', 'thdgrice')
gb = field(glob.glob(RUNS['GLOBAL (11Xg)'])[0], 'thdgrice', 'thdgrice')
if ga is not None and gb is not None:
    n = min(ga.shape[0], gb.shape[0])
    ma, mb = np.nanmean(ga[:n], axis=1), np.nanmean(gb[:n], axis=1)
    print('   GSMART: ' + ' '.join(f'{nm}{v:+8.2e}' for nm, v in zip(names, ma)))
    print('   GLOBAL: ' + ' '.join(f'{nm}{v:+8.2e}' for nm, v in zip(names, mb)))
    print(f'   annual mean  GSMART {np.nanmean(ma):+.4e}   GLOBAL {np.nanmean(mb):+.4e}   '
          f'ratio {np.nanmean(mb)/np.nanmean(ma) if np.nanmean(ma) else float("nan"):.3f}')
    print(f'   extreme |thdgrice| per month, GSMART max {np.nanmax(np.abs(ga)):.3e}  '
          f'GLOBAL max {np.nanmax(np.abs(gb)):.3e}')
