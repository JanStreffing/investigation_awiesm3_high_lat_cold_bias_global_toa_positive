"""Winter sea ice in the two hemispheres: what keeps the Southern Ocean pack loose and thin?

Two global levers in ice_thermo_cpl.F90 would act mostly on Southern Ocean winter ice if
the physics differs between hemispheres the way it is expected to:
  gamma_t  ocean-to-ice heat transfer, Qocnice = (sst - Tf(sss)) * gamma_t * cc
           (gamma_t = 10 m/day, cc = 1025*4190), large only where the water under the
           ice is well above freezing;
  h0       lead closing, dA = dh_openwater / clamp(h, 0.5, 1.5), which matters only if
           the loose pack is thermodynamic rather than opened by the drift.
This measures the inputs to that argument, per hemisphere in local winter (SH JJA,
NH DJF of the same model year), over ice (a_ice >= 0.15), the model pack (>= 0.8) and,
in the south, the cells HadISST has at >= 0.8 but the model below 0.8 ("loose"):
  * concentration tendencies from thermodynamics (thdgrarea) and from the dynamics step
    (dyngrarea: advection, divergence and the cut-off at A = 1, taken across the
    ice_setup_step dynamics call), in frac/month;
  * volume tendencies thdgrice and dyngrice, as cm of grid-mean ice per month;
  * thermal forcing sst - Tf and the implied Qocnice, from daily means;
  * mixed-layer depth MLD2, ice speed (mean of daily |u_ice|), floe thickness m_ice/a_ice.
Plus the SH monthly area budget, integrated over 50-90S, in M km2/month.

FESOM regular 0.5-degree output (*.fesom.gr.*), cos-lat cell areas.  Daily fields:
a_ice sst sss uice vice.  Monthly: m_ice thdgrarea dyngrarea thdgrice dyngrice MLD2.

Usage:  ARM=PI200 Y0=1530 Y1=1539 python3 scripts/analysis/ice_winter_budget_hemispheres.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); Y0, Y1 = int(os.environ.get('Y0', 1530)), int(os.environ.get('Y1', 1539))
SIC = os.environ.get('SIC', '/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc')
GAMMA_T, CC = 10.0 / 86400.0, 1025.0 * 4190.0
MON = 86400.0 * 30.44
NY = Y1 - Y0 + 1


def rd(v, y):
    with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
        a = np.squeeze(d[v].values)
        lat = np.squeeze(d['lat'].values); lon = np.squeeze(d['lon'].values)
    a = np.where(np.abs(a) > 1e30, np.nan, a).astype('f8')
    return a, lat, lon


def monthly(a):
    if a.shape[0] == 12:
        return a
    n = [31, 29 if a.shape[0] == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    e = np.cumsum([0] + n)
    return np.stack([np.nanmean(a[e[i]:e[i + 1]], 0) for i in range(12)])


F = {}
for y in range(Y0, Y1 + 1):
    cur = {}
    for v in ('a_ice', 'm_ice', 'thdgrarea', 'dyngrarea', 'thdgrice', 'dyngrice', 'MLD2'):
        cur[v], lat, lon = rd(v, y)
        cur[v] = monthly(cur[v])
    sst, _, _ = rd('sst', y); sss, _, _ = rd('sss', y)
    cur['dT'] = monthly(sst - (-0.0575 * sss + 1.7105e-3 * np.sqrt(np.maximum(sss, 0) ** 3) - 2.155e-4 * sss ** 2))
    del sst, sss
    u, _, _ = rd('uice', y); w, _, _ = rd('vice', y)
    cur['spd'] = monthly(np.hypot(u, w)); del u, w
    for k, a in cur.items():
        F[k] = F.get(k, 0) + a / NY

LAT = np.broadcast_to(lat[:, None], F['a_ice'].shape[1:])
AREA = np.cos(np.deg2rad(LAT)) * (0.5 * 111.195e3) ** 2
with xr.open_dataset(SIC, decode_times=False) as ds:
    sic = np.squeeze(ds['sic'].values); slat = np.squeeze(ds['latitude'].values); slon = np.squeeze(ds['longitude'].values)
if np.nanmax(sic) > 1.5: sic = sic / 100.0
iy = np.abs(lat[:, None] - slat[None, :]).argmin(1)
ix = np.abs(((lon[:, None] % 360 - slon[None, :] + 180) % 360) - 180).argmin(1)
sic = np.nan_to_num(sic[:, iy, :][:, :, ix], nan=0.0)


def season(ms):
    return {k: np.nanmean(a[ms], 0) for k, a in F.items()}, sic[ms].mean(0)


def row(name, S, k):
    k = k & np.isfinite(S['a_ice']) & np.isfinite(S['dT'])
    if not k.any():
        return
    w = AREA[k]
    av = lambda f: float(np.nansum((f[k] * w)) / np.sum(w[np.isfinite(f[k])]))
    print(f'  {name:<16}{np.sum(w) / 1e12:6.2f}{av(S["a_ice"]):6.2f}'
          f'{av(S["thdgrarea"]) * MON:+9.3f}{av(S["dyngrarea"]) * MON:+9.3f}'
          f'{av(S["thdgrice"]) * MON * 100:+8.1f}{av(S["dyngrice"]) * MON * 100:+8.1f}'
          f'{av(S["dT"]):8.3f}{av(S["dT"]) * GAMMA_T * CC:8.1f}{av(S["MLD2"]):7.0f}'
          f'{av(S["spd"]) * 100:7.1f}{av(S["m_ice"] / np.maximum(S["a_ice"], 0.15)):7.2f}')


print(f'{ARM} {Y0}-{Y1}, local winter.  Area tendencies frac/month, volume tendencies cm/month (grid mean),')
print('dT = sst - Tf [K], Qoi = dT*gamma_t*cc [W/m2], MLD2 [m], ice speed [cm/s], floe thickness [m]')
hdr = f'  {"class":<16}{"Mkm2":>6}{"A":>6}{"dA thm":>9}{"dA dyn":>9}{"dV thm":>8}{"dV dyn":>8}{"dT":>8}{"Qoi":>8}{"MLD":>7}{"spd":>7}{"h":>7}'
for hemi, ms, band in (('SH JJA', [5, 6, 7], LAT < -50), ('NH DJF', [11, 0, 1], LAT > 50)):
    S, obs = season(ms)
    print(f'\n{hemi}'); print(hdr)
    row('ice >= 0.15', S, band & (S['a_ice'] >= 0.15))
    row('pack >= 0.8', S, band & (S['a_ice'] >= 0.8))
    row('0.15-0.8', S, band & (S['a_ice'] >= 0.15) & (S['a_ice'] < 0.8))
    if hemi.startswith('SH'):
        row('loose vs HadISST', S, band & (obs >= 0.8) & (S['a_ice'] >= 0.15) & (S['a_ice'] < 0.8))

print('\nSH 50-90S monthly area budget [M km2/month]: thermodynamic, dynamic, and model/HadISST area')
for m in range(2, 10):
    k = LAT < -50
    t = np.nansum((F['thdgrarea'][m] * AREA)[k]) * MON / 1e12
    d = np.nansum((F['dyngrarea'][m] * AREA)[k]) * MON / 1e12
    a = np.nansum((F['a_ice'][m] * AREA)[k]) / 1e12
    o = np.nansum((sic[m] * AREA)[k & np.isfinite(F['a_ice'][m])]) / 1e12
    print(f'  month {m + 1:2d}  thm {t:+6.2f}  dyn {d:+6.2f}  area model {a:5.2f}  HadISST {o:5.2f}')
