"""Effective ocean-to-ice exchange velocity under McPhee St*u*, against the old constant.

u* is recomputed exactly as ice_thermo_cpl.F90 does it, from DAILY ice and ocean surface
velocity: u* = sqrt(Cd_oce_ice) * |u_ice - u_oc|, Cd_oce_ice = 0.0055 (namelist.ice).
The McPhee arm then exchanges at St* * max(u*, u*_min), St* = 0.0057, u*_min = 0.005 m/s,
against gamma_t = 10 m/day before.  Reported per hemisphere and month over the pack
(a_ice >= 0.15), daily values weighted by daily ice area, with the ratio to 10 m/day and
the actual McPhee flux St* u* (sst - Tf) rho cp computed daily, so the u* dT correlation
is kept.

Usage:  ROOT=<runtime> ARM=PICAL_crunveg_ihf1 Y0=2127 Y1=2129 python ice_ocean_ustar_by_hemisphere.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R, ARM = os.environ['ROOT'], os.environ['ARM']
Y0, Y1 = int(os.environ['Y0']), int(os.environ['Y1'])
CD, ST, UMIN = 0.0055, 0.0057, 0.005
CC = 1025.0 * 4190.0
MN = 'JFMAMJJASOND'


def rd(v, y):
    with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
        a = np.squeeze(d[v].values).astype('f4')
        lat = np.squeeze(d['lat'].values)
    return np.where(np.abs(a) > 1e30, np.nan, a), lat


def month_index(n):
    ml = [31, 29 if n == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return np.repeat(np.arange(12), ml)[:n]


S = {}
for y in range(Y0, Y1 + 1):
    ai, lat = rd('a_ice', y)
    ui, _ = rd('uice', y); vi, _ = rd('vice', y)
    uo, _ = rd('unod_sfc', y); vo, _ = rd('vnod_sfc', y)
    sst, _ = rd('sst', y); sss, _ = rd('sss', y)
    n = ai.shape[0]
    shapes = {k: a.shape[0] for k, a in (('uice', ui), ('vice', vi), ('unod', uo), ('vnod', vo), ('sst', sst), ('sss', sss))}
    assert all(v == n for v in shapes.values()), f'mixed frequencies: a_ice {n}, {shapes}'
    ustar = np.sqrt(CD) * np.hypot(ui - uo, vi - vo)
    dT = sst - (-0.0575 * sss + 1.7105e-3 * np.sqrt(np.maximum(sss, 0) ** 3) - 2.155e-4 * sss ** 2)
    # CONST=1 scores the control arm: exchange at the constant 10 m/day instead
    gam = np.full_like(ustar, 10.0 / 86400.0) if os.environ.get('CONST') == '1' else ST * np.maximum(ustar, UMIN)
    fields = {'us': ustar, 'gam': gam, 'q': gam * dT * CC, 'dT': dT}
    mi = month_index(n)
    for m in range(12):
        k = mi == m
        w = np.nan_to_num(ai[k])
        S.setdefault('ai', np.zeros((12,) + ai.shape[1:]))[m] += np.nanmean(ai[k], 0)
        S.setdefault('w', np.zeros((12,) + ai.shape[1:]))[m] += w.sum(0)
        for name, a in fields.items():
            S.setdefault(name, np.zeros((12,) + ai.shape[1:]))[m] += np.nansum(np.nan_to_num(a[k]) * w, 0)
S['ai'] /= (Y1 - Y0 + 1)

LAT = np.broadcast_to(lat[:, None], S['ai'].shape[1:])
AREA = np.cos(np.deg2rad(LAT))
print(__doc__.split('Usage:')[0])
print(f'{ARM} {Y0}-{Y1}' + ('  (CONST=1: exchange at 10 m/day)' if os.environ.get('CONST') == '1' else '') + '\n')
for hemi, sel in (('NH', LAT > 40), ('NH>80N', LAT > 80), ('SH', LAT < -40), ('SH<65S', LAT < -65)):
    print(f'{hemi:7s} mon   u*[mm/s]  St*u*[m/day]  ratio to 10  Qocnice[W/m2 ice]         sst-Tf[K]')
    for m in range(12):
        k = sel & (S['ai'][m] >= 0.15)
        W = (S['w'][m] * AREA)[k].sum()
        if W == 0:
            continue
        f = lambda nm: (S[nm][m] * AREA)[k].sum() / W
        g = f('gam') * 86400
        print(f'        {MN[m]}   {1e3 * f("us"):8.2f}  {g:12.2f}  {g / 10:11.2f}  {f("q"):23.1f}  {f("dT"):9.3f}')
    print()
