"""Southern Ocean winter ice and snow thickness, PI200 against GIOMAS.

OpenIFS receives FESOM's ice and snow thickness (A_Ice_thickness, A_Snow_thickness in
namcouple) and conducts heat through them on its ice tile, so thin ice or thin snow in
winter is a direct route to a warm skin over the pack.  GIOMAS (1989-2018 climatology) is
a model-data assimilation product, not an observation, and gives only effective ice
thickness.  Compared: mean EFFECTIVE thickness (volume per cell area) over cells with
effective thickness >= 0.05 m, 55-80S, JJA and Sep; for the model also floe thickness
(m_ice/a_ice) and snow on ice (m_snow/a_ice), both where a_ice >= 0.15.

Usage:  ARM=PI200 Y0=1520 Y1=1539 python3 scripts/analysis/so_winter_ice_thickness.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); Y0, Y1 = int(os.environ.get('Y0', 1520)), int(os.environ.get('Y1', 1539))
GIO = '/work/ab0246/a270092/obs/GIOMAS/GIOMAS_heff_miss_time_mon.nc'
SEAS = {'JJA': [5, 6, 7], 'Sep': [8]}


def gr(v, y):
    with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
        a = np.squeeze(d[v].values).astype('f8')
        if a.shape[0] > 12:                                  # a_ice is daily
            n = [31, 29 if a.shape[0] == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            e = np.cumsum([0] + n); a = np.stack([a[e[i]:e[i + 1]].mean(0) for i in range(12)])
        return a, np.squeeze(d['lat'].values)


acc = {}
for y in range(Y0, Y1 + 1):
    for v in ('a_ice', 'm_ice', 'm_snow'):
        a, lat = gr(v, y); acc[v] = acc.get(v, 0) + a / (Y1 - Y0 + 1)
LAT = np.broadcast_to(lat[:, None], acc['a_ice'].shape[1:]); W = np.cos(np.deg2rad(LAT))
band = (LAT >= -80) & (LAT <= -55)
with xr.open_dataset(GIO, decode_times=False) as d:
    k = [c for c in d.data_vars if c.lower().startswith('heff')][0]
    H = np.squeeze(d[k].values).astype('f8')
    glat = np.squeeze(d['lat_scaler'].values); gW = np.squeeze(d['dxt'].values * d['dyt'].values)
gband = (glat >= -80) & (glat <= -55)
print(f'{ARM} {Y0}-{Y1} against GIOMAS 1989-2018, 55-80S')
print(f'  {"":<5}{"eff. thick model":>17}{"GIOMAS":>8}{"floe thick":>11}{"snow on ice":>12}{"mean a_ice":>11}')
for s, m in SEAS.items():
    a = np.nanmean(acc['a_ice'][m], 0); h = np.nanmean(acc['m_ice'][m], 0); sn = np.nanmean(acc['m_snow'][m], 0)
    g = np.nanmean(H[m], 0)
    ke = band & np.isfinite(h) & (h >= 0.05); kg = gband & np.isfinite(g) & (g >= 0.05); ki = band & np.isfinite(a) & (a >= 0.15)
    av = lambda f, k, w=W: float(np.average(f[k], weights=w[k]))
    print(f'  {s:<5}{av(h, ke):17.2f}{av(g, kg, gW):8.2f}{av(h / a, ki):11.2f}{av(sn / a, ki):12.2f}{av(a, ki):11.2f}')
