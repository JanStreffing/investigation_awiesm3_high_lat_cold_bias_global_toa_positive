"""First decade (2171-2179) of the arms queued on 2026-10-06/07 against PICAL_crunveg_tke_albsn082.

Arms, each one change from the control at the 2169 restart:
  icb      Antarctic calving as icebergs (calving - basal melt), mapper calving off
  wamz0    wave model on the ECWAM default roughness (LLGCBZ0, LLNORMAGAM)
  cdoi10   ice-ocean drag 0.0055 -> 0.010
  cpres15  ice strength c_pressure 20 -> 15
Reports, where the output exists: Southern Ocean cap (September density step 200-300 m minus 0-50 m,
surface salinity, 50-150 m temperature) in the open ocean south of 55S deeper than 2000 m and in the
regions NEW / OLD of so_destabilisation.py, mean of 2175-2179; September and February sea-ice area
south of 50S and March / September north of 40N; mapper runoff south of 50S and iceberg melt (Gt/yr);
zonal-mean 10 m wind, eastward stress and evaporation over the ocean at 50-65S; global T2m and net TOA.

Usage: [ARMS=icb,ctl2,... ARMS_TAG=_x] overnight_arms_first_decade.py [y0 y1]
(reval environment; writes data/clim/overnight_arms<ARMS_TAG>_<y0>-<y1>.txt)
"""
import os
import numpy as np
import xarray as xr
import gsw

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
B = f'{P}/runtime/awiesm3-v3.4'
RUNS = [('control', 'PICAL_crunveg_tke_albsn082')] + [(k, f'PICAL_crunveg_{k}') for k in
        os.environ.get('ARMS', 'icb,wamz0,cdoi10,cpres15').split(',')]          # ARMS=icb,ctl2,z0def,... selects the runs
TAG = os.environ.get('ARMS_TAG', '')
import sys
YA, YB = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (2171, 2179)   # atmosphere and ice means
YC = range(max(YA, YB - 4), YB + 1)                                                    # cap: last five Septembers
md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lat = md['lat'].values; area = md['nod_area'].values[0]
zi = np.abs(md['nz'].values); zm, dz = 0.5 * (zi[1:] + zi[:-1]), np.diff(zi)
cav = md['ulevels_nod2D'].values > 1
od = (~cav) & (np.abs(md['zbar_n_bottom'].values) > 2000) & (lat < -55)
mld = lambda e: np.abs(xr.open_dataset(f'{B}/{e}/outdata/fesom/MLD2.fesom.2177.nc')['MLD2'].isel(time=8).values)
ma, mt = mld('PICAL_crunveg_kpp_gm1000'), mld('PICAL_crunveg_tke_albsn082')
REG = {'SO55': od, 'NEW': od & (ma > 1000) & (mt < 500), 'OLD': od & (ma > 1000) & (mt > 1000)}
GT = 1000.0 * 86400 * 365.25 / 1e12


def have(exp, var, y):
    return os.path.isfile(f'{B}/{exp}/outdata/fesom/{var}.fesom.{y}.nc')


def lay(x, a, b):
    m = (zm >= a) & (zm < b)
    return np.nansum(x[:, m] * dz[m], 1) / np.nansum(np.where(np.isfinite(x[:, m]), dz[m], 0), 1)


def cap(exp):
    out = {}
    ys = [y for y in YC if have(exp, 'temp', y) and have(exp, 'salt', y)]
    if not ys:
        return None
    for k, reg in REG.items():
        idx = np.where(reg)[0]; w = area[idx] / area[idx].sum(); acc = []
        for y in ys:
            T = xr.open_dataset(f'{B}/{exp}/outdata/fesom/temp.fesom.{y}.nc')['temp'].isel(time=8).values[idx]
            S = xr.open_dataset(f'{B}/{exp}/outdata/fesom/salt.fesom.{y}.nc')['salt'].isel(time=8).values[idx]
            T = np.where(S > 1, T, np.nan); S = np.where(S > 1, S, np.nan)
            f = lambda x, a, b: float(np.nansum(lay(x, a, b) * w))
            acc.append((float(gsw.sigma0(f(S, 200, 300), f(T, 200, 300)) - gsw.sigma0(f(S, 0, 50), f(T, 0, 50))), f(S, 0, 50), f(T, 50, 150)))
        out[k] = np.mean(acc, 0)
    return out, len(ys)


def ice(exp):
    r = {}
    for nm, sel, mo in (('SH Sep', lat < -50, 8), ('SH Feb', lat < -50, 1), ('NH Mar', lat > 40, 2), ('NH Sep', lat > 40, 8)):
        v = []
        for y in range(YA, YB + 1):
            if not have(exp, 'a_ice', y):
                continue
            a = xr.open_dataset(f'{B}/{exp}/outdata/fesom/a_ice.fesom.{y}.nc')['a_ice']
            a = a.sel(time=a['time'].dt.month == mo + 1).mean('time').values
            v.append(np.nansum((a * area)[sel]) / 1e12)
        r[nm] = np.mean(v) if v else np.nan
    return r


def fw(exp):
    ro, ib, n = 0.0, 0.0, 0
    for y in range(YA, YB + 1):
        if not have(exp, 'runoff', y):
            continue
        r = np.nan_to_num(xr.open_dataset(f'{B}/{exp}/outdata/fesom/runoff.fesom.{y}.nc')['runoff'].mean('time').values)
        ro += np.sum((r * area)[lat < -50]) * GT; n += 1
        for v in ('ibfwb', 'ibfwbv', 'ibfwl', 'ibfwe'):
            if have(exp, v, y):
                x = np.nan_to_num(xr.open_dataset(f'{B}/{exp}/outdata/fesom/{v}.fesom.{y}.nc')[v].mean('time').values)
                ib += np.sum(x * area) * GT
    return (ro / n, ib / n) if n else (np.nan, np.nan)


def atm(exp):
    d = f'{B}/{exp}/outdata/oifs'
    ys = [y for y in range(YA, YB + 1) if os.path.isfile(f'{d}/atm_remapped_1m_ewss_{y}-{y}.nc')]
    if not ys:
        return None

    def tm(v):
        a = 0
        for y in ys:
            x = xr.open_dataset(f'{d}/atm_remapped_1m_{v}_{y}-{y}.nc', decode_times=False)[v]
            a = a + x.mean([k for k in x.dims if k not in ('lat', 'lon')]).load()
        return (a / len(ys)).sortby('lat')
    l = xr.open_dataset(f'{d}/atm_remapped_1m_lsm_{ys[-1]}-{ys[-1]}.nc', decode_times=False)['lsm']
    oc = l.mean([k for k in l.dims if k not in ('lat', 'lon')]).sortby('lat') < 0.5
    u, s, e, t2 = tm('10u'), tm('ewss') / 3600, tm('e') * 24000, tm('2t')
    toa = (tm('tsr') + tm('ttr')) / 3600
    w = np.cos(np.deg2rad(u['lat']))
    band = lambda x, a, b, m=True: float((x * w).where(oc & (x['lat'] >= a) & (x['lat'] <= b)).sum() / (w * xr.ones_like(x)).where(oc & (x['lat'] >= a) & (x['lat'] <= b)).sum())
    gm = lambda x: float((x.mean('lon') * w).sum() / w.sum())
    return dict(u=band(u, -65, -50), tau=band(s, -65, -50), cd=1e3 * band(s, -65, -50) / band(u, -65, -50) ** 2, evap=-band(e, -65, -50),
                utrop=band(u, -20, 20), t2m=gm(t2) - 273.15, toa=gm(toa), n=len(ys))


out = [f'Arms from the 2169 restart against PICAL_crunveg_tke_albsn082 (control). Ice, fluxes, atmosphere: {YA}-{YB}; cap: last five Septembers.',
       'PHC3 September: NEW step 0.320, S 0-50 34.04, T 50-150 -1.12; OLD (Weddell) step 0.238, S 34.22, T -1.10.', '']
C = {k: cap(e) for k, e in RUNS}
out.append(f"{'cap, September':16s}" + ''.join(f'{k:>24s}' for k in REG) + '   years')
out.append(f"{'':16s}" + ''.join(f"{'step   S0-50  T50-150':>24s}" for k in REG))
for k, e in RUNS:
    if C[k] is None:
        out.append(f'{k:16s}   no 2175-2179 output yet'); continue
    out.append(f'{k:16s}' + ''.join(f'{C[k][0][r][0]:10.3f}{C[k][0][r][1]:8.3f}{C[k][0][r][2]:6.2f}' for r in REG) + f'   {C[k][1]}')
out += ['', f"{'sea-ice area 1e6 km2':22s}{'SH Sep':>8s}{'SH Feb':>8s}{'NH Mar':>8s}{'NH Sep':>8s}{'runoff S of 50S':>17s}{'iceberg melt':>14s}  [Gt/yr]"]
for k, e in RUNS:
    i = ice(e); f = fw(e)
    out.append(f"{k:22s}{i['SH Sep']:8.2f}{i['SH Feb']:8.2f}{i['NH Mar']:8.2f}{i['NH Sep']:8.2f}{f[0]:17.0f}{f[1]:14.0f}")
out += ['', f"{'atmosphere':16s}{'u10 50-65S':>11s}{'tau 50-65S':>11s}{'tau/u2 e3':>10s}{'evap 50-65S':>12s}{'u10 20S-20N':>12s}{'T2m glob':>9s}{'net TOA':>8s}{'yrs':>4s}   (ERA5: u10 5.51, tau 0.145, tau/u2 4.78)"]
for k, e in RUNS:
    a = atm(e)
    if a is None:
        out.append(f'{k:16s}   no output yet'); continue
    out.append(f"{k:16s}{a['u']:11.2f}{a['tau']:11.3f}{a['cd']:10.2f}{a['evap']:12.2f}{a['utrop']:12.2f}{a['t2m']:9.2f}{a['toa']:8.2f}{a['n']:4d}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/overnight_arms{TAG}_{YA}-{YB}.txt', 'w').write(txt + '\n')
