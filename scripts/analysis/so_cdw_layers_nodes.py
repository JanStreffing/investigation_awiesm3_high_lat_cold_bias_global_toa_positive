"""60-78S layer temperatures against PHC3 from FESOM node output (runs without regular-grid files).

Same layers and weighting as so_cdw_drift.py (volume-weighted, open-ocean nodes, no cavity),
for node output temp.fesom.YYYY.nc with dims (time, nz1|nz, nod2) or (time, nod2, nz),
annual means every STEP years.  MESH must be the run's mesh diag (nod_area, lat) and PHC the
PHC3 file interpolated to that mesh.  Also prints the mesh spacing at 60-78S.

Usage:  ROOT=/work/bb1469/a270092/runtime/awicm3-v3.3.0 ARM=PI Y0=1850 Y1=2014 STEP=15 \
        MESH=/work/ab0246/a270092/input/fesom2/core2/fesom.mesh.diag.nc \
        PHC=/work/ab0246/a270092/postprocessing/climatologies/CORE2/temp.fesom.1958.nc \
        python3 scripts/analysis/so_cdw_layers_nodes.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ['ROOT']; ARM = os.environ['ARM']; MESH = os.environ['MESH']; PHC = os.environ['PHC']
Y0, Y1, STEP = int(os.environ.get('Y0', 1850)), int(os.environ.get('Y1', 2014)), int(os.environ.get('STEP', 15))
LAYERS = [('0-100', 0, 100), ('100-500', 100, 500), ('500-1500', 500, 1500), ('1500-bot', 1500, 1e5)]
with xr.open_dataset(MESH) as m:
    lat = m['lat'].values; area = m['nod_area'].values[0]; nlv = m['nlevels_nod2D'].values; z = m['nz1'].values
    ul = m['ulevels_nod2D'].values if 'ulevels_nod2D' in m else np.ones(len(lat), int)
NZ = len(z); k = (lat >= -78) & (lat <= -60) & (ul == 1)
edges = np.concatenate([[0.0], 0.5 * (z[1:] + z[:-1]), [z[-1] + 0.5 * (z[-1] - z[-2])]]); dz = np.diff(edges)
print(f'{ARM}: mesh spacing 60-78S {np.sqrt(np.average(area[k], weights=area[k])) / 1e3:.1f} km (area-weighted), {len(lat)} nodes')


def to_nodes_levels(a):
    a = np.squeeze(a).astype('f8')
    if a.ndim == 3: a = a.mean(0)
    if a.shape[0] != len(lat): a = a.T
    a[np.abs(a) > 100] = np.nan
    return np.where(np.arange(a.shape[1])[None, :] < (nlv - 1)[:, None], a[:, :NZ], np.nan)


def layers(a):
    out = {}
    for n, lo, hi in LAYERS:
        kz = (z >= lo) & (z < hi); w = area[k][:, None] * dz[None, kz] * np.isfinite(a[k][:, kz])
        out[n] = np.nansum(np.nan_to_num(a[k][:, kz]) * w) / w.sum()
    return out


with xr.open_dataset(PHC, decode_times=False) as d:
    ref = layers(to_nodes_levels(d['temp'].values[..., :NZ] if d['temp'].shape[-1] > NZ else d['temp'].values))
print('  PHC3: ' + '  '.join(f'{n} {ref[n]:.2f}' for n in ref))
print(f'  {"year":<6}' + ''.join(f'{n:>10}' for n in ref) + '     (model minus PHC3)')
for y in range(Y0, Y1 + 1, STEP):
    f = f'{R}/{ARM}/outdata/fesom/temp.fesom.{y}.nc'
    if not os.path.exists(f): continue
    with xr.open_dataset(f, decode_times=False) as d:
        t = layers(to_nodes_levels(d['temp'].values))
    print(f'  {y:<6}' + ''.join(f'{t[n]:7.2f} ({t[n] - ref[n]:+.2f})' for n in ref))
