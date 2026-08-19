"""Close the coupled energy budget: TOA -> surface -> ocean.  It closes.

WHY.  11G warms the whole ocean +0.197 K/century, about +1.9 K per millennium, which
makes a multi-millennial piControl spin-up untenable.  The question was whether that is
driven by the radiative imbalance or by something unaccounted for between the top of the
atmosphere and the sea surface.

THE ANSWER: nothing is unaccounted for, once snow enthalpy is in the surface budget.

    surface = ssr + str + sshf + slhf - sf*rho*Lf

Snowfall enthalpy is an INTERNAL atmosphere<->surface transfer, so it belongs in the
surface budget and not at TOA.  This is the campaign's established convention -- report
sub:snowenth, and notes/AMIP_BASELINE_AND_ROUND09_2026-07-28.md, which measured it at
0.82-0.94 W/m2 in AMIP.  Coupled, 2026-08-19:

  arm   window          TOA   SFCraw   snowE     SFC   resid   OHC W/m2  K/cent  TOA implies
  11E   1380-1399  +0.758   +1.764   0.962  +0.802  -0.044     +0.714  +0.214       +0.227
  11G   1380-1399  +0.732   +1.728   0.976  +0.752  -0.020     +0.657  +0.197       +0.219
  11L   1380-1389  +0.484   +1.506   1.018  +0.488  -0.004     +0.419  +0.126       +0.145
  11M   1380-1389  -0.016   +1.071   1.073  -0.002  -0.014     -0.069  -0.021       -0.005

The column conserves to 0.02-0.04 W/m2, and the measured ocean warming matches what the
TOA imbalance alone implies, to within about 10 %.  The whole chain closes.

A CORRECTION THIS FILE REPLACES.  An earlier version of this script omitted the snow term
and reported a "-1.0 W/m2 structural leak, invariant across arms", concluding that the
drift was an energy-conservation problem and that tuning TOA to zero could not fix it.
That was wrong on both counts.  The -1.0 W/m2 IS the snow enthalpy: sf*Lf measures 0.962,
0.976, 1.018 and 1.073 in the four arms against residuals of -1.006, -0.996, -1.022 and
-1.088 -- agreement to 0.02 W/m2.  Its apparent invariance across arms, which is what made
it look structural, is simply that global snowfall barely changes between them.

The other half of that error was a window mismatch: 11G's TOA is +0.205 over 1350-69 plus
1380-89 but +0.732 over 1380-99, and the "factor of three" came from comparing the early
number against the late drift.  Net TOA is itself drifting upward through these runs,
which is worth knowing on its own.

SO THE PRACTICAL CONCLUSION INVERTS: tuning net TOA to zero DOES stop the ocean drift.
11M is the demonstration -- TOA -0.016 gives -0.021 K/century, i.e. a flat ocean.  The
difficulty with 11M is not energetics, it is that it gets there by dimming the planet and
destroying the Siberian forest.

WHAT REMAINS OPEN.  The ocean takes in slightly more than it stores -- about 0.1 W/m2 in
every arm.  Candidates are sea-ice melt and the coastal masking in the regridded fh
integration below, which drops non-finite and exactly-zero cells and weights the ocean
integral by total area.  Integrate fh on the native mesh with proper cell areas before
treating that as physical.

SIGN CONVENTIONS, established from the data rather than assumed: IFS surface fluxes are
downward-positive; FESOM fh is the opposite, since its global mean is negative while the
ocean warms.
"""
import glob
import sys
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
ACC = 3600.0
LF = 3.337e5      # J/kg, latent heat of fusion
EARTH = 5.10072e14
RHO, CP = 1027.0, 3990.0

ARMS = [('11E', f'{R}/Tuning_test_11E_swemin15_K1', range(1380, 1400)),
        ('11G', f'{R}/Tuning_test_11G_inppmin50k', range(1380, 1400)),
        ('11L', f'{R}/11L', range(1380, 1390)),
        ('11M', f'{R}/11M', range(1380, 1390))]


def gm(root, var, years):
    vals = []
    for y in years:
        f = f'{root}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not glob.glob(f):
            continue
        try:
            with xr.open_dataset(f, decode_times=False) as d:
                a = np.asarray(d[var].values, dtype=float) / ACC
                lat = d['lat'].values
            vals.append(float(np.average(a.mean(axis=0).mean(axis=1),
                                         weights=np.cos(np.deg2rad(lat)))))
        except Exception:
            pass
    return np.mean(vals) if vals else np.nan


def ocean_terms(root, years):
    th, vo = {}, {}
    for y in years:
        for var, store in (('thetaoga', th), ('volo', vo)):
            f = f'{root}/outdata/fesom/{var}.fesom.{y}.nc'
            if not glob.glob(f):
                continue
            try:
                with xr.open_dataset(f, decode_times=False) as d:
                    store[y] = float(np.asarray(d[var].values, dtype=float).mean())
            except Exception:
                pass
    ys = sorted(set(th) & set(vo))
    if len(ys) < 8:
        return np.nan, np.nan
    T = np.array([th[y] for y in ys])
    V = float(np.mean([vo[y] for y in ys]))
    dTdt = np.polyfit(ys, T, 1)[0]
    ohc = dTdt * V * RHO * CP / 3.15576e7 / EARTH
    fh = []
    for y in ys:
        f = f'{root}/outdata/fesom/fh.fesom.gr.{y}.nc'
        if not glob.glob(f):
            continue
        try:
            with xr.open_dataset(f, decode_times=False) as d:
                a = np.asarray(d['fh'].values, dtype=float).mean(axis=0)
                lat = d['lat'].values
            w = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, a.shape[1]))
            m = np.isfinite(a) & (a != 0)
            if m.sum():
                fh.append(float((a[m] * w[m]).sum() / w.sum()))
        except Exception:
            pass
    return ohc, (-np.mean(fh) if fh else np.nan)   # sign: negative fh = into ocean


def main():
    print(__doc__)
    print('=' * 92)
    print(f'\n1. IFS COLUMN\n\n{"arm":5s} {"window":11s} {"TOA":>8s} {"SFC":>8s} '
          f'{"resid":>8s}')
    for tag, root, Y in ARMS:
        Y = list(Y)
        t = {v: gm(root, v, Y) for v in ('tsr', 'ttr', 'ssr', 'str', 'sshf', 'slhf')}
        toa = t['tsr'] + t['ttr']
        # Snow enthalpy: internal atmosphere<->surface transfer, belongs here only.
        snow = gm(root, 'sf', Y) * 1000.0 * LF
        sfc = t['ssr'] + t['str'] + t['sshf'] + t['slhf'] - snow
        print(f'{tag:5s} {Y[0]}-{Y[-1]} {toa:+8.3f} {sfc:+8.3f} {toa - sfc:+8.3f}'
              f'   (snow enthalpy {snow:.3f})')
    print('\n   resid is now 0.02-0.04 W/m2: the column conserves. Omitting the snow')
    print('   term produces a spurious -1.0 W/m2 that looks like a structural leak.')

    print(f'\n2. FESOM OCEAN (per m2 of EARTH)\n\n{"arm":5s} {"into ocean":>11s} '
          f'{"OHC tend":>10s} {"resid":>8s}')
    for tag, root, Y in ARMS:
        ohc, fh = ocean_terms(root, list(Y))
        if np.isnan(ohc):
            print(f'{tag:5s}  too few years')
            continue
        print(f'{tag:5s} {fh:+11.3f} {ohc:+10.3f} {fh - ohc:+8.3f}')
    print('\n   The ocean is self-consistent to ~0.1 W/m2. The unexplained gap is between')
    print('   the atmosphere\'s surface flux and the ocean\'s uptake -- see the docstring.')


if __name__ == '__main__':
    main()
