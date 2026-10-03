"""Does a Fox-Kemper mixed-layer-eddy parameterisation make sense on CORE3?

CORE3 is variable resolution, so "our resolution" is not one number.  The test for whether
an MLE parameterisation is warranted is whether the model resolves the instability it
represents: baroclinic instability of mixed-layer fronts, whose scale is the MIXED-LAYER
deformation radius

    Ld_ML = N_ML * H / |f|          N_ML = rms buoyancy frequency over the mixed layer

If Ld_ML is unresolved, the model cannot generate the restratification itself and a
parameterisation is doing real work.  For comparison the FIRST BAROCLINIC radius

    Ld_1  = c_1 / |f|               c_1 = (1/pi) * integral(N dz)

is what GM is scaled against, and FESOM already carries `mesh_resolution(n)`, so
per-node resolution scaling on this mesh is existing machinery rather than new work.

The Fox-Kemper coefficient itself carries an explicit resolution factor, Delta_s/L_f,
because the grid-scale buoyancy gradient underestimates the true frontal one.  On a
variable mesh that factor is a FIELD, not a constant, which is the implementation point
this script is meant to size.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

D = '/work/bb1469/a270092/runtime/awiesm3-v3.4/11X/outdata/fesom'
MESHF = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
YEARS = range(1380, 1390)
OMEGA = 7.2921e-5
L_F = 4000.0          # frontal width assumed by Fox-Kemper et al. (2011), metres

BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45), ('45-30S', -45, -30),
         ('30S-30N', -30, 30), ('30-45N', 30, 45), ('45-60N', 45, 60), ('60-90N', 60, 90)]
BANDS = BANDS[::-1]          # report north to south

m = xr.open_dataset(MESHF, decode_times=False)
lat = m['lat'].values.astype('f8')
area = m['cell_area'].values.astype('f8')
assert lat.size == 220509, 'MESH GUARD: not core3'
# effective grid spacing: FESOM's mesh_resolution is ~sqrt of the scalar cell area
res = np.sqrt(area)
f = 2 * OMEGA * np.sin(np.deg2rad(lat))
absf = np.maximum(np.abs(f), 1e-5)          # regularised near the equator

with xr.open_dataset(f'{D}/temp.fesom.1389.nc', decode_times=False) as d:
    zmid = d['nz'].values.astype('f8')
dbn = m['depth_bnds'].values.astype('f8')
nz = zmid.size
dz = np.abs(dbn[1:nz + 1]) - np.abs(dbn[:nz])


def dec(var, months=None):
    acc, n = None, 0
    for y in YEARS:
        with xr.open_dataset(f'{D}/{var}.fesom.{y}.nc', decode_times=False) as d:
            v = [k for k in d.data_vars if k not in
                 ('bounds_lon', 'bounds_lat', 'time_bounds')][-1]
            a = d[v].values.astype('f8')
        a = a[months].mean(axis=0) if months is not None else a.mean(axis=0)
        acc = a if acc is None else acc + a
        n += 1
    return acc / n


# local winter: NH Feb-Mar, SH Aug-Sep -- the deepest mixed layers
n2_nh = dec('N2', [1, 2]); n2_sh = dec('N2', [7, 8])
ml_nh = np.abs(dec('MLD2', [1, 2])); ml_sh = np.abs(dec('MLD2', [7, 8]))
north = lat >= 0
N2 = np.where(north[:, None], n2_nh[:, :nz], n2_sh[:, :nz])
H = np.where(north, ml_nh, ml_sh)

# N over the mixed layer, and over the full column for the first baroclinic mode
N = np.sqrt(np.maximum(N2, 1e-10))
inml = (zmid[None, :] <= H[:, None]) & np.isfinite(N2)
w = np.where(inml, dz[None, :], 0.0)
N_ml = np.where(w.sum(axis=1) > 0,
                np.nansum(N * w, axis=1) / np.maximum(w.sum(axis=1), 1e-30), np.nan)
wet = np.isfinite(N2)
c1 = np.nansum(N * dz[None, :] * wet, axis=1) / np.pi

Ld_ml = N_ml * H / absf
Ld_1 = c1 / absf
ok = np.isfinite(Ld_ml) & (H > 10)

print('=' * 104)
print('CORE3 RESOLUTION vs THE SCALES A MIXED-LAYER-EDDY SCHEME WOULD REPRESENT')
print('11X 1380-89, local winter (NH Feb-Mar, SH Aug-Sep); area-weighted band means')
print('=' * 104)
print(f'{"band":<11}{"grid dx":>10}{"(min-max)":>17}{"winter H":>10}{"Ld_ML":>9}'
      f'{"dx/Ld_ML":>10}{"Ld_1":>9}{"dx/Ld_1":>9}{"dx/L_f":>8}')
print(f'{"":11}{"[km]":>10}{"[km]":>17}{"[m]":>10}{"[km]":>9}{"":10}{"[km]":>9}{"":9}{"":8}')
for nm, lo, hi in BANDS:
    s = ok & (lat >= lo) & (lat < hi)
    if not s.any():
        continue
    aw = lambda x: np.sum(x[s] * area[s]) / np.sum(area[s])
    r = aw(res) / 1e3
    lm = aw(Ld_ml) / 1e3
    l1 = aw(Ld_1) / 1e3
    print(f'{nm:<11}{r:>10.1f}{f"{res[s].min()/1e3:.0f}-{res[s].max()/1e3:.0f}":>17}'
          f'{aw(H):>10.0f}{lm:>9.1f}{r/max(lm,1e-9):>10.1f}{l1:>9.1f}'
          f'{r/max(l1,1e-9):>9.1f}{aw(res)/L_F:>8.1f}')
s = ok
aw = lambda x: np.sum(x[s] * area[s]) / np.sum(area[s])
print(f'{"GLOBAL":<11}{aw(res)/1e3:>10.1f}{f"{res[s].min()/1e3:.0f}-{res[s].max()/1e3:.0f}":>17}'
      f'{aw(H):>10.0f}{aw(Ld_ml)/1e3:>9.1f}{aw(res)/aw(Ld_ml):>10.1f}'
      f'{aw(Ld_1)/1e3:>9.1f}{aw(res)/aw(Ld_1):>9.1f}{aw(res)/L_F:>8.1f}')

print('\n  dx/Ld_ML  >> 1 means mixed-layer baroclinic instability is unresolved,')
print('            i.e. the model cannot make this restratification itself.')
print('  dx/L_f    is the Fox-Kemper resolution factor (L_f = 4 km frontal width).')

# how much does the resolution factor vary across the mesh?
print('\n' + '-' * 104)
print('THE VARIABLE-MESH PROBLEM: how far off would a single tuned coefficient be?')
print('-' * 104)
q = np.percentile(res[ok] / 1e3, [1, 5, 25, 50, 75, 95, 99])
print('  grid dx percentiles [km]: ' + '  '.join(f'p{p}={v:.0f}' for p, v in
      zip([1, 5, 25, 50, 75, 95, 99], q)))
print(f'  ratio p99/p1 = {q[-1]/q[0]:.1f}  -> a constant Delta_s/L_f would be wrong by that '
      f'factor across the mesh')
outc = ok & (((lat >= 30) & (lat < 60)) | ((lat >= -60) & (lat < -30)))
print(f'  at the outcrops (30-60N and 30-60S): dx = {np.sum(res[outc]*area[outc])/np.sum(area[outc])/1e3:.1f} km '
      f'(range {res[outc].min()/1e3:.0f}-{res[outc].max()/1e3:.0f})')
print(f'  refscalresol used by GM is 100 km, so GM there is scaled by '
      f'~{(np.sum(res[outc]*area[outc])/np.sum(area[outc])/1e5):.2f}')
