"""Model minus PHC3 ocean temperature by latitude band and depth: what constrains the uptake?

WHY.  The heat uptake is stored at 700-2000 m and the polar bands are the doorway it enters
through.  Before reaching for a mixing lever to change how deep polar uptake penetrates, the
question is whether the model's polar water column is already wrong, and in which direction
at which depth.  That is an observational question and PHC3 is the reference the model was
INITIALISED from, so any departure is drift the model produced itself, not an initialisation
offset.  This is the cleanest constraint available here.

REFERENCE.  /work/ab0246/a270092/postprocessing/climatologies/CORE3/temp.fesom.1958.nc, the
PHC3 climatology already interpolated to the CORE3 mesh, 211567 nodes, the same mesh as the
model output.  NOT obs/phc3/temp.fesom.1958.nc, which has 126858 nodes: that file is on a
different mesh and cannot be differenced against this model at all.

WHAT PHC3 IS AND IS NOT.  A hydrographic climatology built mostly from pre-Argo bottle and
CTD data, with the Arctic filled from Russian sources.  It is a MEAN STATE, not a time
series, so it constrains where the water column sits, not how fast heat is taken up.  Polar
coverage below 1000 m and under ice is thin and the deep Southern Ocean is the weakest part
of it.  Treat the deep polar rows as indicative.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESH = '/work/ab0246/a270092/input/fesom2/core3_beta'
REF = '/work/ab0246/a270092/postprocessing/climatologies/CORE3/temp.fesom.1958.nc'

BANDS = [('90-60S  Antarctic/SO', -90, -60),
         ('60-45S  subantarctic', -60, -45),
         ('30S-30N  tropics', -30, 30),
         ('45-60N  subpolar NA', 45, 60),
         ('60-90N  Arctic/Nordic', 60, 90)]
DEPTHS = [('0-100 m', 0, 100), ('100-700 m', 100, 700),
          ('700-2000 m', 700, 2000), ('>2000 m', 2000, 1e9)]

ARM = os.environ.get('ARM', '11Q')
YEARS = [int(y) for y in os.environ.get('YEARS', '1382,1383,1384,1385,1386,1387,1388,1389').split(',')]

root = f'{R092}/{ARM}' if os.path.isdir(f'{R092}/{ARM}') else \
       [p for p in glob.glob(f'{R092}/*{ARM}*') if os.path.isdir(p)][0]

with xr.open_dataset(f'{MESH}/fesom.mesh.diag.nc', decode_times=False) as m:
    lat = m['lat'].values
    if np.abs(lat).max() < 4:
        lat = np.rad2deg(lat)
    nod_area = m['nod_area'].values
    zbar = m['nz'].values
z = np.abs(zbar); dz = np.diff(z); nlev = len(dz)
zmid = 0.5 * (z[:-1] + z[1:])
vol = nod_area[:nlev, :] * dz[:, None]

with xr.open_dataset(REF, decode_times=False) as d:
    n = 'temp' if 'temp' in d.data_vars else list(d.data_vars)[-1]
    ref = np.squeeze(d[n].values)
if ref.shape[0] != nlev:
    ref = ref.T
ref = ref[:nlev, :]
if ref.shape[1] != len(lat):
    raise SystemExit(f'  mesh mismatch: reference has {ref.shape[1]} nodes, mesh has {len(lat)}')

acc, n_ok = None, 0
for y in YEARS:
    f = glob.glob(f'{root}/outdata/fesom/temp.fesom.{y}.nc')
    if not f:
        continue
    with xr.open_dataset(f[0], decode_times=False) as d:
        nm = 'temp' if 'temp' in d.data_vars else list(d.data_vars)[-1]
        a = d[nm].mean(dim='time').values.T.astype(np.float64)[:nlev, :]
    acc = a if acc is None else acc + a
    n_ok += 1
if not n_ok:
    raise SystemExit('  no model temp output found')
mod = acc / n_ok
print(f'arm {ARM}, {n_ok} yr ({YEARS[0]}-{YEARS[-1]}), reference PHC3 on CORE3\n')

good = np.isfinite(mod) & np.isfinite(ref) & (vol > 0)
d3 = np.where(good, mod - ref, np.nan)
w3 = np.where(good, vol, 0.0)

hdr = f'  {"band":24s}' + ''.join(f'{d[0]:>12s}' for d in DEPTHS)
print('  MODEL minus PHC3, volume-weighted mean temperature [K]\n')
print(hdr); print('  ' + '-' * (len(hdr) - 2))
for name, la0, la1 in BANDS:
    sel = (lat >= la0) & (lat < la1)
    row = []
    for _dn, dd0, dd1 in DEPTHS:
        lk = (zmid >= dd0) & (zmid < dd1)
        ww = w3[np.ix_(lk, sel)]; vv = d3[np.ix_(lk, sel)]
        row.append(np.nansum(vv * ww) / ww.sum() if ww.sum() else np.nan)
    print(f'  {name:24s}' + ''.join(f'{v:+12.3f}' for v in row))
row = []
for _dn, dd0, dd1 in DEPTHS:
    lk = (zmid >= dd0) & (zmid < dd1)
    ww = w3[lk, :]; vv = d3[lk, :]
    row.append(np.nansum(vv * ww) / ww.sum() if ww.sum() else np.nan)
print('  ' + '-' * (len(hdr) - 2))
print(f'  {"GLOBAL":24s}' + ''.join(f'{v:+12.3f}' for v in row))
