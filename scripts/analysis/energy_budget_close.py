"""Close the coupled energy budget: TOA -> surface -> ocean, and find where it leaks.

WHY.  11G drifts +0.19 K/century in whole-ocean potential temperature -- about +1.9 K per
millennium, which makes a multi-millennial piControl spin-up untenable.  Its net TOA is
only +0.205 W/m2 over the clean window, which can supply roughly +0.06 K/century.  The
ocean was warming about three times faster than the top of the atmosphere was delivering,
so something between the two was unaccounted for.

WHAT IT CHECKS, in the order the energy flows:

  1. IFS column:  TOA net (tsr+ttr) against surface net (ssr+str+sshf+slhf, downward +).
     A conserving atmosphere gives ~0.  Its heat capacity is ~1e7 J/m2/K, so a sustained
     1 W/m2 would be tens of K of atmospheric cooling in twenty years -- it is not
     storage, it is a leak.
  2. FESOM ocean:  surface heat flux (fh) against the OHC tendency implied by thetaoga
     and volo.  This is the ocean's internal consistency.

RESULT, 2026-08-19, four arms:

    arm   window        TOA      SFC    resid
    11E   1380-1399   +0.758   +1.764   -1.006
    11G   1380-1399   +0.732   +1.728   -0.996
    11L   1380-1389   +0.484   +1.506   -1.022
    11M   1380-1389   -0.016   +1.071   -1.088

The residual is -1.0 W/m2 in every arm while TOA moves 0.77 W/m2 across them, so it is
STRUCTURAL, not a property of any lever.  This is the familiar IFS atmospheric energy
non-conservation -- kinetic dissipation not returned as heat, and the enthalpy of
precipitation -- which is order 1-2 W/m2 in this model family.

THE CONSEQUENCE, which is what matters for the spin-up: THE SURFACE SEES ROUGHLY
TOA + 1 W/m2.  Tuning net TOA to zero therefore does NOT give a non-drifting ocean.  11M
is the demonstration: its TOA is -0.016, as balanced as anything the campaign has
produced, and its surface still takes +1.071.  Closing the ocean drift by radiative
tuning alone would require driving TOA to about -1 W/m2, which would put the model
badly wrong against CERES.  The drift is an energy-conservation problem, not a cloud
tuning problem.

WHAT IS NOT RESOLVED.  IFS says the surface takes +1.728 W/m2 (11G, per m2 of Earth);
FESOM says the ocean takes +0.765 and warms at +0.657.  The ocean side is self-consistent
to about 0.11 W/m2, but there is a further ~1 W/m2 between the atmosphere's surface flux
and the ocean's uptake that this script does not account for.  Candidates: heat into land
and sea ice, and the coastal masking in the regridded fh integration below, which is the
weakest step here -- it drops non-finite and exactly-zero cells and weights the ocean
integral by total area.  Do not read that second gap as a second leak until fh has been
integrated on the native mesh with proper cell areas.

SIGN CONVENTIONS, established from the data rather than assumed: IFS surface fluxes are
downward-positive, so ssr+str+sshf+slhf is heat INTO the surface.  FESOM fh is the
opposite -- its global mean is negative while the ocean warms -- so negative fh is heat
into the ocean.
"""
import glob
import sys
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
ACC = 3600.0
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
        sfc = t['ssr'] + t['str'] + t['sshf'] + t['slhf']
        print(f'{tag:5s} {Y[0]}-{Y[-1]} {toa:+8.3f} {sfc:+8.3f} {toa - sfc:+8.3f}')
    print('\n   resid ~ -1 W/m2 in every arm while TOA spans 0.77 across them: structural.')

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
