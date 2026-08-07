"""Scale ladder: how must the depletion curve change with model resolution?

WHY THIS EXISTS.  Mode 3 has NO scale awareness.  d_c, SCALE and SWEMIN are pure
TCO95 calibrations, and mode 2's SDOR term was never carried into mode 3.  A subgrid
depletion curve represents variance the grid CANNOT resolve, so as the grid refines
the parameterised depletion must weaken toward a step function -- applying TCO95
constants at TCO4000 (~2.5 km) would double-count patchiness the model now resolves
explicitly.  AWI-ESM3 is expected to run from TCO95 up to about TCO4000, a factor ~40
in grid spacing, so this is not a corner case.

WHAT IS ALREADY KNOWN, at the two ENDS of that range:
  * RIHMI snow courses are 1-2 km transects -- essentially TCO4000 resolution.  The
    fitted curve (SCALE=1, SWEMIN at the bare numerical floor of 3) IS the
    high-resolution limit, measured directly.
  * Rutgers aggregated to the TCO95 Siberian box says SWEMIN ~ 30 at ~100 km.
Two measured anchors, but nothing in between and no functional form.

WHY SWEMIN IS THE RIGHT CARRIER.  It sets the minimum box-mean snow MASS before cover
may saturate.  At 2 km a thin uniform dusting genuinely does cover everything, so the
floor should be negligible.  At 100 km a box-mean SWE of 3 kg/m2 is a few snowy cells
among many bare ones, so far more box-mean mass is needed before the whole box fills.
"How much box-mean mass before the box is full" is inherently a box-SIZE quantity.
That it also happens to be the numerical stability floor is an accident.

METHOD, and why it is paired rather than climatological.
  amip_presentday carries monthly sd/rsn for 1989-2014, which OVERLAPS Rutgers.  So
  for every model cell and month we can pair the model's own box-mean depth and
  density with the observed covered fraction from the SAME period -- an empirical
  SCF(d, rho) at TCO95, not a comparison of two climatologies from different eras.
  Rutgers 24 km binary cells are area-averaged into each model cell.

  The ladder then coarse-grains BOTH sides into n x n blocks of model cells, giving
  the empirical curve at ~1x, 2x, 3x, 4x the TCO95 spacing, and SWEMIN is refitted at
  each scale with every other parameter held at the P3 values.  Because both sides
  are coarsened identically, the only thing that changes is the box size.

LIMITS, stated because they bound the extrapolation:
  * Rutgers' own cells are 24 km aggregates, so nothing here constrains BELOW ~24 km.
    The 2.5 km end still rests entirely on the snow courses.
  * The ladder only goes COARSER than TCO95; the refinement direction is inferred by
    extrapolating toward the course anchor, not measured.
  * Siberian patchiness is mostly non-orographic (wind, vegetation, aspect), which is
    exactly why mode 2's SDOR term could carry only ~16 % of the scale.  A single
    power law may therefore not hold outside this domain.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

RUT = '/work/ab0246/a270092/obs/snowcover/rutgers_nh_24km_weekly_sce.nc'
RUN = 'amip_presentday'
YEARS = range(1990, 2015)
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
BASE = dict(dcl=0.014, dch=0.026, md=4.70, bl=1.46, bh=0.40,
            dcmax=0.30, rhoref=200.0, scale=1.0)
# months where the curve does any work; Jul/Aug are snow-free, Dec-Mar saturated
FITMON = [8, 9, 10, 3, 4, 5]          # Sep Oct Nov Apr May Jun


def scf(depth, rho, cvh, swemin):
    p = BASE
    r = np.maximum(rho, 50.0) / p['rhoref']
    floor = swemin / np.maximum(rho, 1.0)
    dcl = np.minimum(p['dcmax'], np.maximum(p['scale'] * p['dcl'] * r ** p['md'], floor))
    dch = np.minimum(p['dcmax'], np.maximum(p['scale'] * p['dch'] * r ** p['md'], floor))
    sl = np.clip((depth / np.maximum(dcl, 1e-9)) ** p['bl'], 0, 1)
    sh = np.clip((depth / np.maximum(dch, 1e-9)) ** p['bh'], 0, 1)
    f = np.clip(cvh, 0.0, 1.0)
    live = (depth > 1e-6) & (depth * rho > 1e-6)
    return np.where(live, np.clip((1 - f) * sl + f * sh, 0, 1), 0.0)


print(__doc__.split('LIMITS,')[0])
print('=' * 98)

# ---- model grid and Siberian box ---------------------------------------------
with xr.open_dataset(LSMF) as d:
    lsm = d['lsm'].isel(time_counter=0).values
    mlat, mlon = d['lat'].values, d['lon'].values
dlat = abs(float(mlat[1] - mlat[0])); dlon = abs(float(mlon[1] - mlon[0]))
km_per_deg = 111.32
# The remapped grid is ANISOTROPIC in km: 0.938 deg meridional is ~104 km while
# 0.900 deg zonal is only ~42 km at 65N.  The box SIZE that matters for a subgrid
# fraction is the equivalent-area length sqrt(dy*dx), not the zonal spacing -- using
# the zonal number alone understated the scale by 1.6x in the first run.
dy_km = dlat * km_per_deg
dxz_km = dlon * km_per_deg * np.cos(np.deg2rad(65.0))
dx_km = float(np.sqrt(dy_km * dxz_km))
print(f'model grid: {dlat:.3f} x {dlon:.3f} deg -> {dy_km:.0f} km merid x {dxz_km:.0f} km '
      f'zonal at 65N; equivalent-area length {dx_km:.0f} km')
print('NOTE: SCF is computed on the NATIVE TCO95 grid (~100 km), not on this remapped\n'
      '      output grid, so the n=1 rung is FINER than the model actually works at.\n')

iy0 = int(np.argmin(np.abs(mlat - 75))); iy1 = int(np.argmin(np.abs(mlat - 55)))
if iy0 > iy1:
    iy0, iy1 = iy1, iy0
ix0 = int(np.argmin(np.abs(mlon - 60))); ix1 = int(np.argmin(np.abs(mlon - 180)))
SY, SX = slice(iy0, iy1 + 1), slice(ix0, ix1 + 1)
sub_lat, sub_lon = mlat[SY], mlon[SX]
land = lsm[SY, SX] >= 0.5
print(f'box: {land.shape[0]} x {land.shape[1]} model cells, {land.sum()} on land')

# ---- Rutgers -> model cells, monthly climatology 1990-2014 -------------------
ds = xr.open_dataset(RUT, decode_times=True)
rlat = ds['latitude'].values; rlon = ds['longitude'].values
rland = ds['land'].values; rarea = ds['area'].values
ok = np.isfinite(rlat) & (np.abs(rlat) <= 90) & (rland > 0)
inbox = ok & (rlat >= sub_lat.min()) & (rlat <= sub_lat.max()) & \
        (((rlon + 360) % 360) >= 60) & (((rlon + 360) % 360) <= 180)
ry, rx = np.where(inbox)
# nearest model cell for each Rutgers cell
jy = np.abs(sub_lat[None, :] - rlat[ry, rx][:, None]).argmin(axis=1)
jx = np.abs(sub_lon[None, :] - ((rlon[ry, rx] + 360) % 360)[:, None]).argmin(axis=1)
flat_idx = jy * land.shape[1] + jx
wgt = rarea[ry, rx].astype('float64')
ncell = land.size
print(f'Rutgers cells mapped into the box: {ry.size}')

t = ds['time'].values
yr = np.array([int(str(x)[:4]) for x in t])
mo = np.array([int(str(x)[5:7]) for x in t])
sel_t = np.flatnonzero((yr >= 1990) & (yr <= 2014))

obs_num = np.zeros((12, ncell)); obs_den = np.zeros((12, ncell))
for i in sel_t:
    sce = ds['snow_cover_extent'].isel(time=i).values[ry, rx].astype('float64')
    good = np.isfinite(sce) & (sce <= 1)
    m = mo[i] - 1
    np.add.at(obs_num[m], flat_idx[good], (sce * wgt)[good])
    np.add.at(obs_den[m], flat_idx[good], wgt[good])
ds.close()
OBS = np.where(obs_den > 0, obs_num / np.maximum(obs_den, 1e-12), np.nan).reshape(12, *land.shape)
print(f'observed cover built from {sel_t.size} weekly fields\n')

# ---- model monthly depth / density, same period ------------------------------
def monthly(var):
    acc, n = None, 0
    for y in YEARS:
        f = f'{RT}/{RUN}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values[:, SY, SX]
        acc = a if acc is None else acc + a
        n += 1
    return acc / n, n


sd, n1 = monthly('sd')
rsn, n2 = monthly('rsn')
rho = np.maximum(rsn, 1e-6)
depth = sd * 1000.0 / rho
print(f'model depth/density from {RUN}, {n1} yr')

# high-vegetation fraction is static in AMIP; take it from a P3 daily file
cvh = np.zeros_like(depth[0])
for y in (1900,):
    f = f'{RT}/amip_P3_scffit/outdata/oifs/atm_remapped_1d_cvh_1d_{y}-{y}.nc'
    if os.path.exists(f):
        with xr.open_dataset(f, decode_times=False) as d:
            cvh = d['cvh'].values[0][SY, SX]
print(f'cvh: mean {np.nanmean(cvh[land]):.3f} over box land\n')

# ---- the ladder ---------------------------------------------------------------
def block(a, n, w):
    """Area-weighted n x n block mean of a (…, ny, nx) field; partial blocks dropped."""
    ny, nx = a.shape[-2:]
    ny2, nx2 = (ny // n) * n, (nx // n) * n
    a = a[..., :ny2, :nx2]; w = w[:ny2, :nx2]
    sh = a.shape[:-2] + (ny2 // n, n, nx2 // n, n)
    aw = (a * w).reshape(sh).sum(axis=(-3, -1))
    ww = w.reshape(ny2 // n, n, nx2 // n, n).sum(axis=(-3, -1))
    return aw / np.maximum(ww, 1e-12), ww


coslat = np.cos(np.deg2rad(np.broadcast_to(sub_lat[:, None], land.shape)))
wland = np.where(land, coslat, 0.0)

print('SWEMIN refitted at each box size, all other parameters held at the P3 values\n')
print(f'  {"blocks":>7s}{"~dx km":>9s}{"ncells":>8s}{"SWEMIN":>9s}{"RMSE":>9s}'
      f'{"RMSE@3":>9s}{"improve":>9s}')
cand = np.arange(0.5, 400.1, 0.5)
ladder = []
for n in (1, 2, 3, 4, 6):
    D, W = block(depth, n, wland)
    R, _ = block(rho, n, wland)
    C, _ = block(cvh, n, wland)
    O, _ = block(np.where(np.isfinite(OBS), OBS, np.nan), n, wland)
    m = np.isfinite(O) & (W[None, :, :] > 0)
    keep = np.zeros_like(m); keep[FITMON] = True
    m &= keep
    if m.sum() < 50:
        continue
    d_, r_, c_, o_ = D[m], R[m], np.broadcast_to(C[None], D.shape)[m], O[m]
    ww = np.broadcast_to(W[None], D.shape)[m]
    best = min(((np.sqrt(np.sum(ww * (scf(d_, r_, c_, s) - o_) ** 2) / ww.sum()), s)
                for s in cand))
    r3 = np.sqrt(np.sum(ww * (scf(d_, r_, c_, 3.0) - o_) ** 2) / ww.sum())
    ladder.append((n, n * dx_km, best[1], best[0]))
    print(f'  {n:7d}{n*dx_km:9.0f}{m.sum():8d}{best[1]:9.1f}{best[0]:9.4f}'
          f'{r3:9.4f}{100*(r3-best[0])/r3:8.0f}%')

print(f'\n  snow-course anchor (RIHMI, 1-2 km transects): SWEMIN ~ 3 (the numerical floor)')
if len(ladder) >= 2:
    L = np.array([x[1] for x in ladder]); S = np.array([x[2] for x in ladder])
    good = S > 0
    if good.sum() >= 2:
        p = np.polyfit(np.log(L[good]), np.log(S[good]), 1)
        print(f'  power-law fit over the measured ladder: SWEMIN ~ {np.exp(p[1]):.3g} '
              f'* dx_km^{p[0]:.2f}')
        for target, name in ((dx_km, 'TCO95'), (25.0, 'TCO400'), (2.5, 'TCO4000')):
            print(f'    -> {name:8s} (dx {target:6.1f} km): SWEMIN = '
                  f'{np.exp(p[1]) * target ** p[0]:6.1f}')

print("""

  READING IT.  If SWEMIN rises monotonically with box size and the fitted exponent is
  positive, the scale hypothesis holds and SWEMIN is the carrier.  Compare the
  extrapolation to TCO4000 against the snow-course anchor of ~3: if they agree, the
  power law spans the whole range AWI-ESM3 needs and can go into the code as a
  function of grid spacing.  If the extrapolation undershoots or overshoots badly,
  the law is only valid coarser than TCO95 and the fine end must stay a documented
  namelist choice rather than an automatic formula.""")
