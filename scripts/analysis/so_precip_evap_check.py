"""Atmospheric freshwater input to the Southern Ocean: precipitation and evaporation against GPCP.

Zonal mean over the ocean of total precipitation (cp + lsp) and evaporation (e) of
PICAL_crunveg_tke_albsn082, 2195-2209, by latitude band, with GPCP v2.3 precipitation (1979-2021 mean,
obs/gpcp/precip.mon.mean.nc) on the model grid under the model's ocean mask. P - E is the largest
freshwater source of the surface cap south of 55S; evaporation has no observational counterpart here.
GPCP is poorly constrained over sea ice and at high southern latitudes (few gauges, satellite
retrievals over ice), so differences of 10-20 % there are within its uncertainty.

Usage: so_precip_evap_check.py   (reval environment; writes data/clim/so_precip_evap_check.txt)
"""
import os
import numpy as np
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082/outdata/oifs'
Y0, Y1 = 2195, 2209
MMD = 1000.0 * 24.0                 # m of water per hour -> mm/day


def tm(v):
    a = 0
    for y in range(Y0, Y1 + 1):
        d = xr.open_dataset(f'{R}/atm_remapped_1m_{v}_{y}-{y}.nc', decode_times=False)[v]
        a = a + d.mean([k for k in d.dims if k not in ('lat', 'lon')]).load()
    return (a / (Y1 - Y0 + 1)).sortby('lat')


P = (tm('cp') + tm('lsp')) * MMD
E = tm('e') * MMD                   # negative = evaporation
l = xr.open_dataset(f'{R}/atm_remapped_1m_lsm_{Y1}-{Y1}.nc', decode_times=False)['lsm']
oc = l.mean([k for k in l.dims if k not in ('lat', 'lon')]).sortby('lat') < 0.5
g = xr.open_dataset('/albedo/work/user/jstreffi/obs/gpcp/precip.mon.mean.nc')['precip'].mean('time').sortby('lat')
G = g.interp(lat=P['lat'], lon=P['lon'], kwargs=dict(fill_value='extrapolate'))
w = np.cos(np.deg2rad(P['lat']))
out = [f'Zonal mean over the ocean [mm/day], PICAL_crunveg_tke_albsn082 {Y0}-{Y1}; GPCP v2.3 1979-2021.',
       f"{'band':10s}{'P model':>9s}{'P GPCP':>9s}{'ratio':>7s}{'E model':>9s}{'P-E model':>11s}{'P-E mm/yr':>11s}"]
for a, b in ((-75, -65), (-65, -55), (-55, -45), (-45, -30), (-30, 30), (30, 60), (60, 90)):
    m = oc & (P['lat'] >= a) & (P['lat'] < b)
    f = lambda x: float((x * w).where(m).sum() / (w * xr.ones_like(x)).where(m).sum())
    p, gp, e = f(P), f(G), f(E)
    out.append(f"{f'{a}..{b}':10s}{p:9.2f}{gp:9.2f}{p / gp:7.2f}{-e:9.2f}{p + e:11.2f}{365.25 * (p + e):11.0f}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/so_precip_evap_check.txt', 'w').write(txt + '\n')
