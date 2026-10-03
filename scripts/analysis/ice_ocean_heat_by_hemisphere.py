"""Would gamma_t (ice-ocean heat transfer) act mostly on the Southern Ocean pack?

In the OpenIFS-coupled build FESOM grows ice from dhice = A*(Qatmice - Qocnice)
(ice_thermo_cpl.F90, coupled-slab convention): the atmosphere owns the surface balance and
the conduction, so FESOM's con/consn do not set growth, but the ocean-to-ice heat flux
    Qocnice = (sst - Tf(sss)) * gamma_t * rho*cp,   gamma_t = 10 m/day (a Fortran parameter)
does, linearly.  A global change of gamma_t is hemispherically selective by physics if the
water under the ice is further above freezing in one hemisphere.  This measures, per
hemisphere and month, over the pack (a_ice >= 0.15): the thermal forcing sst - Tf, the
implied Qocnice (W/m2 of ice area), the ice-area-weighted atmospheric growth it opposes
(thdgrice as growth rate), and the pack area.

FESOM regular 0.5-degree output (*.fesom.gr.*).  Daily or monthly fields are both reduced
to monthly means.

Usage:  ROOT=<runtime> ARM=PICAL_crunveg Y0=2110 Y1=2119 python ice_ocean_heat_by_hemisphere.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = os.environ['ROOT']
ARM = os.environ['ARM']
Y0, Y1 = int(os.environ['Y0']), int(os.environ['Y1'])
GAMMA_T, CC = 10.0 / 86400.0, 1025.0 * 4190.0
MN = 'JFMAMJJASOND'


def monthly(a):
    if a.shape[0] == 12:
        return a
    n = [31, 29 if a.shape[0] == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    e = np.cumsum([0] + n)
    return np.stack([np.nanmean(a[e[i]:e[i + 1]], 0) for i in range(12)])


def rd(v, y):
    with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
        a = np.squeeze(d[v].values).astype('f8')
        lat = np.squeeze(d['lat'].values)
    return monthly(np.where(np.abs(a) > 1e30, np.nan, a)), lat


acc = {}
n = 0
for y in range(Y0, Y1 + 1):
    try:
        ai, lat = rd('a_ice', y)
        sst, _ = rd('sst', y)
        sss, _ = rd('sss', y)
        tg, _ = rd('thdgrice', y)
    except (OSError, KeyError) as e:
        print(f'skip {y}: {e}')
        continue
    tf = -0.0575 * sss + 1.7105e-3 * np.sqrt(np.maximum(sss, 0) ** 3) - 2.155e-4 * sss ** 2
    cur = {'ai': ai, 'dT': sst - tf, 'tg': tg}
    for k, a in cur.items():
        acc[k] = acc.get(k, 0) + np.nan_to_num(a)
    n += 1
F = {k: a / n for k, a in acc.items()}
LAT = np.broadcast_to(lat[:, None], F['ai'].shape[1:])
AREA = np.cos(np.deg2rad(LAT)) * (0.5 * 111.195e3) ** 2

print(__doc__.split('Usage:')[0])
print(f'{ARM} {Y0}-{Y1} ({n} years); pack = a_ice >= 0.15; ice-area weighted\n')
for hemi, sel in (('NH', LAT > 40), ('NH>80N', LAT > 80), ('SH', LAT < -40), ('SH<65S', LAT < -65)):
    print(f'{hemi}  mon  pack[Mkm2]  sst-Tf[K]  Qocnice[W/m2 ice]  thdgrice[cm/mon grid-mean]')
    for m in range(12):
        a = F['ai'][m]
        k = sel & (a >= 0.15)
        if k.sum() == 0:
            continue
        w = (AREA * a)[k]
        dT = np.average(F['dT'][m][k], weights=w)
        q = dT * GAMMA_T * CC
        g = np.average(F['tg'][m][k], weights=AREA[k]) * 86400 * 30.44 * 100
        print(f'    {MN[m]}  {w.sum() / 1e12:9.2f}  {dT:9.3f}  {q:16.1f}  {g:14.2f}')
    print()
