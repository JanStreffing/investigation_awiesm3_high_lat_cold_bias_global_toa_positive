"""Is the warm Circumpolar Deep Water at 60-78S a drift or an equilibrium feature of the run?

Area-weighted annual-mean temperature of 60-78S ocean cells in 0-100, 100-500, 500-1500 and
1500 m-bottom, every STEP years over Y0-Y1, from FESOM's regular output, against PHC3 (mesh
nodes binned to the same bands).  A layer still warming at the end of a 190-year run means
the reservoir under the winter ice grows with time and any ice-side tuning done against it
is done against a moving target.

Usage:  ARM=PI200 Y0=1390 Y1=1579 STEP=5 python3 scripts/analysis/so_cdw_drift.py
        ROOT=/work/bb1469/a270089/runtime/awiesm3-v3.4.2 ARM=AWI-ESM3-VEG-HR-CMIP7-piControl MESH=... PHC=... (HR, DARS2)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); Y0, Y1, STEP = int(os.environ.get('Y0', 1390)), int(os.environ.get('Y1', 1579)), int(os.environ.get('STEP', 5))
LAYERS = [('0-100', 0, 100), ('100-500', 100, 500), ('500-1500', 500, 1500), ('1500-bot', 1500, 1e5)]
MESH = os.environ.get('MESH', '/work/ab0246/a270092/input/fesom2/core3/fesom.mesh.diag.nc')
PHC = os.environ.get('PHC', '/work/ab0246/a270092/postprocessing/climatologies/CORE3_220509/temp.fesom.1958.nc')
# PHC3 on the DARS2 mesh: MESH=.../dars2/fesom.mesh.diag.nc PHC=.../climatologies/DARS2/temp.fesom.1958.nc (57 levels, last one junk)
with xr.open_dataset(MESH) as m:
    nlat = m['lat'].values; area = m['nod_area'].values[0]; ul = m['ulevels_nod2D'].values if 'ulevels_nod2D' in m else np.ones(m.dims['nod2'], int); nlv = m['nlevels_nod2D'].values; z = m['nz1'].values
with xr.open_dataset(PHC, decode_times=False) as d:
    NZ = len(z); p = np.squeeze(d['temp'].values)[:, :NZ].astype('f8'); p[np.abs(p) > 1e30] = np.nan
p = np.where(np.arange(NZ)[None, :] < (nlv - 1)[:, None], p, np.nan)
k = (nlat >= -78) & (nlat <= -60) & (ul == 1)
edges = np.concatenate([[0.0], 0.5 * (z[1:] + z[:-1]), [z[-1] + 0.5 * (z[-1] - z[-2])]]); dz = np.diff(edges)
ref = {}
for n, a, b in LAYERS:
    kz = (z >= a) & (z < b); w = area[k][:, None] * dz[None, kz] * np.isfinite(p[k][:, kz])
    ref[n] = np.nansum(p[k][:, kz] * w) / w.sum()
print(f'{ARM}, 60-78S open ocean, annual mean T [degC] by layer; PHC3: ' + '  '.join(f'{n} {ref[n]:.2f}' for n in ref))
print(f'  {"year":<6}' + ''.join(f'{n:>10}' for n in ref) + '     (model minus PHC3)')
for y in range(Y0, Y1 + 1, STEP):
    f = f'{R}/{ARM}/outdata/fesom/temp.fesom.gr.{y}.nc'
    if not os.path.exists(f): continue
    with xr.open_dataset(f, decode_times=False) as d:
        lat = d['lat'].values; j = np.where((lat >= -78) & (lat <= -60))[0]
        a = d['temp'].isel(lat=slice(j[0], j[-1] + 1)).mean('time').values.astype('f8'); zz = d['nz'].values; la = lat[j]
    a[np.abs(a) > 100] = np.nan                      # fill values differ between runs (1e30 on CORE3, ~1e22 on DARS2 output)
    ed = np.concatenate([[0.0], 0.5 * (zz[1:] + zz[:-1]), [zz[-1] + 0.5 * (zz[-1] - zz[-2])]]); dzz = np.diff(ed)
    ar = np.cos(np.deg2rad(la))[:, None, None] * np.ones(a.shape[1])[None, :, None]
    row = f'  {y:<6}'
    for n, lo, hi in LAYERS:
        kz = (zz >= lo) & (zz < hi); w = ar * dzz[None, None, kz] * np.isfinite(a[:, :, kz])
        t = np.nansum(a[:, :, kz] * w) / w.sum(); row += f'{t:7.2f} ({t - ref[n]:+.2f})'
    print(row)
