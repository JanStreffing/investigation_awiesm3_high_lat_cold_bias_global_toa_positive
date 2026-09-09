"""Where does the heat go?  40-year drift of ocean temperature by depth and band.

THE QUESTION.  The model runs a positive TOA imbalance while the surface stays cold, so
heat is leaving the surface downward.  The mixed layer depth says the tropics are not the
route (mld_baseline_vs_woa18.py: 30S-30N is 8 m too SHALLOW), and the Southern Ocean's
runaway convection is a heat RELEASE pathway, not a burial one -- convection connects the
abyss to a cold surface and makes the deep too cold, which is what the PHC3 comparison
already showed (model -0.40 C at 500 m against PHC3's +1.14).

So the heat has to be going down somewhere else.  This script does not assume where.  It
measures the drift of the model's own temperature between its first and last decade,
resolved by depth and latitude band, and converts it to W/m2 so it can be set against the
TOA imbalance directly.

WHY DRIFT AND NOT BIAS.  A bias against PHC3 mixes what the model inherited from its
initial condition with what it did afterwards.  The difference between this run's first
and last decade is unambiguously what THIS run put there.  (It also sidesteps the mesh
problem: the only PHC3 on disk is on core3_beta, 211567 nodes, and 11X is core3, 220509.)

HEAT CONTENT.  rho*cp*integral(dT dz dA) over each band and depth range, divided by the
elapsed time and by area, giving W/m2.  Two normalisations are printed: per unit area of
the band itself (how hard that band is being loaded) and spread over the global ocean
(that band's share of the global budget, directly comparable to a TOA number).

Usage:  ARM=11X EARLY=1350-1359 LATE=1380-1389 python3 scripts/analysis/ohc_drift_by_depth_band.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESHF = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
N_NODES_CORE3 = 220509
RHO_CP = 1026.0 * 3996.0        # J m-3 K-1
SEC_PER_YEAR = 365.0 * 86400.0

ARM = os.environ.get('ARM', '11X')
_e = os.environ.get('EARLY', '1350-1359').split('-')
_l = os.environ.get('LATE', '1380-1389').split('-')
EARLY = list(range(int(_e[0]), int(_e[-1]) + 1))
LATE = list(range(int(_l[0]), int(_l[-1]) + 1))

BANDS = [('90-60S  Antarctic/SO', -90, -60), ('60-45S  subantarctic', -60, -45),
         ('45-30S  S subtropical', -45, -30), ('30S-30N  tropics', -30, 30),
         ('30-45N  N subtropical', 30, 45), ('45-60N  subpolar NA', 45, 60),
         ('60-90N  Arctic/Nordic', 60, 90)]
SLABS = [('0-100 m', 0, 100), ('100-300 m', 100, 300), ('300-700 m', 300, 700),
         ('700-2000 m', 700, 2000), ('below 2000 m', 2000, 1e9)]


def decade_mean(root, years, var):
    """Volume-resolved decadal mean field (nod2, nz), NaN below the bottom."""
    acc, n = None, 0
    for y in years:
        f = f'{root}/outdata/fesom/{var}.fesom.{y}.nc'
        if not os.path.exists(f):
            print(f'  {y}: missing {var}')
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values.astype('f8').mean(axis=0)      # annual mean of 12 months
        acc = a if acc is None else acc + a
        n += 1
    if n == 0:
        sys.exit(f'no {var} years found')
    return acc / n, n


def main():
    root = f'{R092}/{ARM}'
    with xr.open_dataset(MESHF, decode_times=False) as m:
        lat = m['lat'].values.astype('f8')
        area = m['cell_area'].values.astype('f8')
        dbnds = m['depth_bnds'].values.astype('f8')
    if lat.size != N_NODES_CORE3:
        sys.exit(f'MESH GUARD: {lat.size} nodes, expected {N_NODES_CORE3}')

    with xr.open_dataset(f'{root}/outdata/fesom/temp.fesom.{LATE[-1]}.nc',
                         decode_times=False) as d:
        zmid = d['nz'].values.astype('f8')
    nz = zmid.size
    ztop = np.abs(dbnds[:nz])
    zbot = np.abs(dbnds[1:nz + 1])
    dz = zbot - ztop
    print(f'{nz} levels, top {ztop[0]:.1f}-{zbot[0]:.1f} m, '
          f'deepest {ztop[-1]:.0f}-{zbot[-1]:.0f} m')

    print(f'{ARM}: early {EARLY[0]}-{EARLY[-1]}, late {LATE[0]}-{LATE[-1]}')
    Ta, na = decade_mean(root, EARLY, 'temp')
    Tb, nb = decade_mean(root, LATE, 'temp')
    Sa, _ = decade_mean(root, EARLY, 'salt')
    Sb, _ = decade_mean(root, LATE, 'salt')
    dT = Tb - Ta
    dS = Sb - Sa
    yrs = 0.5 * (LATE[0] + LATE[-1]) - 0.5 * (EARLY[0] + EARLY[-1])
    print(f'  {na} + {nb} years, {yrs:.0f} yr between decade centres')

    wet = np.isfinite(Ta) & np.isfinite(Tb)                   # (nod2, nz)
    vol = area[:, None] * dz[None, :] * wet                   # m3
    global_ocean_area = area[np.any(wet, axis=1)].sum()
    print(f'  global ocean area {global_ocean_area:.4e} m2, '
          f'volume {vol.sum():.4e} m3')

    # ---------------------------------------------------------------- profiles
    print('\n' + '=' * 104)
    print(f'TEMPERATURE DRIFT [K over {yrs:.0f} yr], volume-weighted band means by depth')
    print('=' * 104)
    hdr = f'{"depth [m]":<14}' + ''.join(f'{n.split()[0]:>11}' for n, _, _ in BANDS) + f'{"GLOBAL":>11}'
    print(hdr)
    sel = {}
    for name, lo, hi in BANDS:
        sel[name] = (lat >= lo) & (lat < hi)
    for k in range(nz):
        row = f'{ztop[k]:>6.0f}-{zbot[k]:<7.0f}'
        vals = []
        for name, _, _ in BANDS:
            s = sel[name][:, None] & wet[:, k:k + 1]
            w = (vol[:, k:k + 1] * s).ravel()
            vals.append(np.sum(dT[:, k] * w) / max(w.sum(), 1e-30))
        w = vol[:, k]
        g = np.sum(dT[:, k] * w) / max(w.sum(), 1e-30)
        if k % 2 == 0 or zbot[k] < 300:
            print(row + ''.join(f'{v:>+11.3f}' for v in vals) + f'{g:>+11.3f}')

    # ------------------------------------------------------------ heat content
    print('\n' + '=' * 104)
    print('HEAT CONTENT CHANGE, per unit area OF THE BAND  [W m-2]')
    print('(how hard that band is being loaded)')
    print('=' * 104)
    print(f'{"band":<24}' + ''.join(f'{n:>14}' for n, _, _ in SLABS) + f'{"total":>10}')
    tot_global = np.zeros(len(SLABS))
    for name, lo, hi in BANDS:
        s = sel[name]
        ba = area[s & np.any(wet, axis=1)].sum()
        vals = []
        for j, (sn, z0, z1) in enumerate(SLABS):
            lev = (ztop >= z0) & (ztop < z1)
            q = np.nansum(dT[s][:, lev] * vol[s][:, lev]) * RHO_CP
            tot_global[j] += q
            vals.append(q / ba / (yrs * SEC_PER_YEAR))
        print(f'{name:<24}' + ''.join(f'{v:>+14.3f}' for v in vals) +
              f'{sum(vals):>+10.3f}')
    print('-' * 104)
    print(f'{"GLOBAL (per global ocean area)":<24}' +
          ''.join(f'{q / global_ocean_area / (yrs * SEC_PER_YEAR):>+14.3f}'
                  for q in tot_global) +
          f'{tot_global.sum() / global_ocean_area / (yrs * SEC_PER_YEAR):>+10.3f}')

    print('\n' + '=' * 104)
    print('SAME HEAT, as each band\'s SHARE of the global ocean budget  [W m-2 of global ocean]')
    print('(these add up to the global total, so they can be set against a TOA number)')
    print('=' * 104)
    print(f'{"band":<24}' + ''.join(f'{n:>14}' for n, _, _ in SLABS) + f'{"total":>10}')
    grand = 0.0
    for name, lo, hi in BANDS:
        s = sel[name]
        vals = []
        for sn, z0, z1 in SLABS:
            lev = (ztop >= z0) & (ztop < z1)
            q = np.nansum(dT[s][:, lev] * vol[s][:, lev]) * RHO_CP
            vals.append(q / global_ocean_area / (yrs * SEC_PER_YEAR))
        grand += sum(vals)
        print(f'{name:<24}' + ''.join(f'{v:>+14.3f}' for v in vals) +
              f'{sum(vals):>+10.3f}')
    print('-' * 104)
    print(f'{"TOTAL":<24}' + ' ' * (14 * len(SLABS)) + f'{grand:>+10.3f}')

    # ------------------------------------------------------------------ salt
    print('\n' + '=' * 104)
    print(f'SALINITY DRIFT [psu over {yrs:.0f} yr] -- the SO convection control')
    print('=' * 104)
    print(f'{"band":<24}' + ''.join(f'{n:>14}' for n, _, _ in SLABS))
    for name, lo, hi in BANDS:
        s = sel[name]
        vals = []
        for sn, z0, z1 in SLABS:
            lev = (ztop >= z0) & (ztop < z1)
            w = vol[s][:, lev]
            vals.append(np.nansum(dS[s][:, lev] * w) / max(np.nansum(w), 1e-30))
        print(f'{name:<24}' + ''.join(f'{v:>+14.4f}' for v in vals))

    np.savez(os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), 'data',
        f'ohc_drift_{ARM}.npz'),
        dT=dT, dS=dS, ztop=ztop, zbot=zbot, lat=lat, area=area, dz=dz, yrs=yrs)
    print(f'\nsaved data/ohc_drift_{ARM}.npz')


if __name__ == '__main__':
    main()
