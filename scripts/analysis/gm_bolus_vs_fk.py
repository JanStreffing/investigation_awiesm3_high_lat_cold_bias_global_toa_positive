"""Is FESOM's GM already doing the mixed-layer restratification a Fox-Kemper scheme would?

WHY THIS IS THE DECIDING TEST.  GM runs inside the mixed layer in this configuration --
the line that would zero it there, oce_fer_gm.F90:454, is commented out -- but what it
actually does there is throttled by neutral-slope clipping (k_gm_ktaper=.false., so the
SLOPE is clipped rather than the coefficient).  Two outcomes, and they point opposite ways:

  GM's mixed-layer streamfunction is near zero  -> the gap is real, an MLE scheme has room,
                                                   and implementing Fox-Kemper is warranted.
  GM's is comparable to what Fox-Kemper would add -> we would mostly be swapping one scheme
                                                   for another, and relaxing the slope
                                                   clipping is the far cheaper experiment.

WHAT IS COMPARED.  Both as streamfunctions in m2/s, at the base of the mixed layer:

  Psi_GM = | integral of the GM bolus velocity from the surface to the mixed layer base |
           (bolus_u/v are model output, so this is what the model actually did, not a
           reconstruction from the coefficient)

  Psi_FK = C_e * (dx/L_f) * H^2 * |grad_b| / sqrt(f^2 + tau^-2)      Fox-Kemper et al. 2011

with grad_b the mixed-layer-averaged buoyancy gradient computed on each element from its
three nodes, exactly the way FESOM computes horizontal gradients, and dx the LOCAL grid
scale -- CORE3 spans 5 to 148 km, so this factor is a field, not a constant.

Local winter only (NH Feb-Mar, SH Aug-Sep): restratification matters when the mixed layer
is deep.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, seawater as sw, warnings
warnings.filterwarnings('ignore')

D = '/work/bb1469/a270092/runtime/awiesm3-v3.4/11X/outdata/fesom'
MESHF = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
YEARS = range(1382, 1390)
R_E = 6371000.0
OMEGA = 7.2921e-5
G = 9.81
RHO0 = 1025.0
C_E = 0.06            # Fox-Kemper et al. (2011) coefficient
L_F = 4000.0          # assumed frontal width [m]
TAU = 5 * 86400.0     # frontal lifetime, regularises f near the equator
NH, SH = [1, 2], [7, 8]

BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45), ('45-30S', -45, -30),
         ('30S-30N', -30, 30), ('30-45N', 30, 45), ('45-60N', 45, 60), ('60-90N', 60, 90)]
BANDS = BANDS[::-1]          # report north to south

m = xr.open_dataset(MESHF, decode_times=False)
lon, lat = m['lon'].values.astype('f8'), m['lat'].values.astype('f8')
tri = m['triag_nodes'].values.astype(int) - 1
dbn = m['depth_bnds'].values.astype('f8')
assert lat.size == 220509, 'MESH GUARD: not core3'

with xr.open_dataset(f'{D}/temp.fesom.1389.nc', decode_times=False) as d:
    zmid = d['nz'].values.astype('f8')
nz = zmid.size
dz = np.abs(dbn[1:nz + 1]) - np.abs(dbn[:nz])

# ---------------------------------------------------- element geometry and gradients
elon, elat = lon[tri], lat[tri]
lon0 = elon[:, 0:1]
dlon = ((elon - lon0 + 180.0) % 360.0) - 180.0          # dateline-safe
lat0 = elat.mean(axis=1, keepdims=True)
X = np.deg2rad(dlon) * R_E * np.cos(np.deg2rad(lat0))
Y = np.deg2rad(elat - lat0) * R_E
twoA = ((X[:, 1] - X[:, 0]) * (Y[:, 2] - Y[:, 0]) -
        (X[:, 2] - X[:, 0]) * (Y[:, 1] - Y[:, 0]))
eArea = 0.5 * np.abs(twoA)
ecLat = elat.mean(axis=1)
dx_e = np.sqrt(2.0 * eArea)                              # element grid scale
f_e = 2 * OMEGA * np.sin(np.deg2rad(ecLat))


def grad_on_elem(fn):
    """Horizontal gradient of a node field, on elements."""
    b = fn[tri]
    gx = (b[:, 0] * (Y[:, 1] - Y[:, 2]) + b[:, 1] * (Y[:, 2] - Y[:, 0]) +
          b[:, 2] * (Y[:, 0] - Y[:, 1])) / twoA
    gy = (b[:, 0] * (X[:, 2] - X[:, 1]) + b[:, 1] * (X[:, 0] - X[:, 2]) +
          b[:, 2] * (X[:, 1] - X[:, 0])) / twoA
    return gx, gy


def winter_node(var, nlev):
    """Decade mean of local-winter months, node field (nod2, nlev) or (nod2,)."""
    accN = accS = None
    n = 0
    for y in YEARS:
        with xr.open_dataset(f'{D}/{var}.fesom.{y}.nc', decode_times=False) as d:
            k = [x for x in d.data_vars if x not in
                 ('bounds_lon', 'bounds_lat', 'time_bounds')][-1]
            aN = d[k].isel(time=NH).values.astype('f8').mean(axis=0)
            aS = d[k].isel(time=SH).values.astype('f8').mean(axis=0)
        accN = aN if accN is None else accN + aN
        accS = aS if accS is None else accS + aS
        n += 1
    accN, accS = accN / n, accS / n
    north = lat >= 0
    if accN.ndim == 1:
        return np.where(north, accN, accS)
    return np.where(north[:, None], accN, accS)


def winter_elem3d(var):
    """Decade mean of local-winter months for an ELEMENT 3-D field."""
    accN = accS = None
    n = 0
    for y in YEARS:
        with xr.open_dataset(f'{D}/{var}.fesom.{y}.nc', decode_times=False) as d:
            aN = d[var].isel(time=NH).values.astype('f8').mean(axis=0)
            aS = d[var].isel(time=SH).values.astype('f8').mean(axis=0)
        accN = aN if accN is None else accN + aN
        accS = aS if accS is None else accS + aS
        n += 1
    accN, accS = accN / n, accS / n
    north = ecLat >= 0
    return np.where(north[:, None], accN, accS)


print('loading winter fields ...', flush=True)
T = winter_node('temp', nz)
S = winter_node('salt', nz)
H = np.abs(winter_node('MLD2', 1))
bu = winter_elem3d('bolus_u')
bv = winter_elem3d('bolus_v')
print('  done', flush=True)

# ------------------------------------------------------ mixed-layer buoyancy gradient
sig = sw.dens0(S, T) - 1000.0
b = -G * sig / RHO0
inml = (zmid[None, :] <= np.maximum(H, zmid[0])[:, None]) & np.isfinite(b)
w = np.where(inml, dz[None, :], 0.0)
b_ml = np.where(w.sum(axis=1) > 0,
                np.nansum(b * w, axis=1) / np.maximum(w.sum(axis=1), 1e-30), np.nan)
ok_node = np.isfinite(b_ml)
b_ml_f = np.where(ok_node, b_ml, 0.0)
gx, gy = grad_on_elem(b_ml_f)
gradb = np.sqrt(gx ** 2 + gy ** 2)

# ROBUSTNESS: |grad b| averaged over triangles picks up grid-scale noise, and Psi_FK is
# linear in it, so the raw band mean can exceed published PEAK values.  Recompute with one
# pass of neighbour averaging on b_ml -- if Psi_FK collapses, the signal was noise.
nnl = m['node_node_links'].values.astype(int) - 1
valid = (nnl >= 0)
bs = np.where(ok_node, b_ml, np.nan)
nb = np.where(valid, bs[np.clip(nnl, 0, lat.size - 1)], np.nan)
b_sm = np.nanmean(np.concatenate([bs[:, None], nb], axis=1), axis=1)
b_sm = np.where(np.isfinite(b_sm), b_sm, 0.0)
gxs, gys = grad_on_elem(b_sm)
gradb_sm = np.sqrt(gxs ** 2 + gys ** 2)
good_e = ok_node[tri].all(axis=1) & np.isfinite(gradb) & (eArea > 0)

H_e = H[tri].mean(axis=1)
denom = np.sqrt(f_e ** 2 + TAU ** -2)
Psi_FK = C_E * (dx_e / L_F) * H_e ** 2 * gradb / denom
Psi_FK_sm = C_E * (dx_e / L_F) * H_e ** 2 * gradb_sm / denom

# --------------------------------------------- GM bolus streamfunction at the ML base
# Psi(z) = integral of the bolus velocity from the surface down to z
cum_u = np.nancumsum(np.where(np.isfinite(bu), bu, 0.0) * dz[None, :], axis=1)
cum_v = np.nancumsum(np.where(np.isfinite(bv), bv, 0.0) * dz[None, :], axis=1)
Psi_GM_z = np.sqrt(cum_u ** 2 + cum_v ** 2)
kml = np.clip(np.searchsorted(zmid, H_e), 0, nz - 1)
Psi_GM_ml = Psi_GM_z[np.arange(len(kml)), kml]
Psi_GM_max = np.nanmax(Psi_GM_z, axis=1)

print('\n' + '=' * 104)
print('GM BOLUS STREAMFUNCTION vs WHAT FOX-KEMPER WOULD ADD  [m2 s-1], local winter')
print(f'11X {YEARS[0]}-{YEARS[-1]}, element-area-weighted band means')
print('=' * 104)
print(f'{"band":<11}{"dx":>8}{"H":>7}{"|grad b|":>12}{"Psi_GM at":>12}{"Psi_GM":>10}'
      f'{"Psi_FK":>10}{"FK/GM":>9}{"GM(ml)/":>10}')
print(f'{"":11}{"[km]":>8}{"[m]":>7}{"[s-2]":>12}{"ML base":>12}{"col max":>10}'
      f'{"":10}{"":9}{"GM(max)":>10}')
for nm, lo, hi in BANDS:
    s = good_e & (ecLat >= lo) & (ecLat < hi)
    if s.sum() < 100:
        continue
    aw = lambda x: np.sum(x[s] * eArea[s]) / np.sum(eArea[s])
    gm_ml, gm_mx, fk = aw(Psi_GM_ml), aw(Psi_GM_max), aw(Psi_FK)
    print(f'{nm:<11}{aw(dx_e)/1e3:>8.1f}{aw(H_e):>7.0f}{aw(gradb):>12.3e}'
          f'{gm_ml:>12.4f}{gm_mx:>10.4f}{fk:>10.4f}'
          f'{fk/max(gm_ml,1e-9):>9.1f}{gm_ml/max(gm_mx,1e-9):>10.2f}')
s = good_e
aw = lambda x: np.sum(x[s] * eArea[s]) / np.sum(eArea[s])
gm_ml, gm_mx, fk = aw(Psi_GM_ml), aw(Psi_GM_max), aw(Psi_FK)
print(f'{"GLOBAL":<11}{aw(dx_e)/1e3:>8.1f}{aw(H_e):>7.0f}{aw(gradb):>12.3e}'
      f'{gm_ml:>12.4f}{gm_mx:>10.4f}{fk:>10.4f}{fk/max(gm_ml,1e-9):>9.1f}'
      f'{gm_ml/max(gm_mx,1e-9):>10.2f}')
print('\n  Psi_GM at ML base is what GM actually transports across the mixed layer base.')
print('  GM(ml)/GM(max) << 1 means GM does its work BELOW the mixed layer, not inside it.')

print('\n' + '=' * 104)
print('ROBUSTNESS: is Psi_FK real, or grid-scale noise in |grad b|?')
print('medians alongside means, and Psi_FK recomputed from a once-smoothed buoyancy field')
print('=' * 104)
print(f'{"band":<11}{"Psi_GM mean":>13}{"median":>9}{"Psi_FK mean":>13}{"median":>9}'
      f'{"FK smoothed":>13}{"ratio smth":>12}{"FK_sm/GM":>10}')
for nm, lo, hi in BANDS:
    s2 = good_e & (ecLat >= lo) & (ecLat < hi)
    if s2.sum() < 100:
        continue
    a = lambda x: np.sum(x[s2] * eArea[s2]) / np.sum(eArea[s2])
    md = lambda x: np.median(x[s2])
    print(f'{nm:<11}{a(Psi_GM_ml):>13.4f}{md(Psi_GM_ml):>9.4f}{a(Psi_FK):>13.4f}'
          f'{md(Psi_FK):>9.4f}{a(Psi_FK_sm):>13.4f}'
          f'{a(Psi_FK_sm)/max(a(Psi_FK),1e-9):>12.2f}'
          f'{a(Psi_FK_sm)/max(a(Psi_GM_ml),1e-9):>10.1f}')
s2 = good_e
a = lambda x: np.sum(x[s2] * eArea[s2]) / np.sum(eArea[s2])
md = lambda x: np.median(x[s2])
print(f'{"GLOBAL":<11}{a(Psi_GM_ml):>13.4f}{md(Psi_GM_ml):>9.4f}{a(Psi_FK):>13.4f}'
      f'{md(Psi_FK):>9.4f}{a(Psi_FK_sm):>13.4f}{a(Psi_FK_sm)/max(a(Psi_FK),1e-9):>12.2f}'
      f'{a(Psi_FK_sm)/max(a(Psi_GM_ml),1e-9):>10.1f}')
