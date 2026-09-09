"""Does FESOM's MLD2 measure the same thing as the WOA18 mixed-layer climatology?

WHY THIS EXISTS.  The campaign wants to know whether the model mixes surface heat down
too efficiently, and the only mixing quantity with a real observational counterpart is
the mixed layer depth.  A model-minus-obs map is worthless if the two numbers are not the
same quantity, so this audits the definitions before any comparison is made.

THE TWO DEFINITIONS, from the source rather than from memory.

  WOA18 (woa18_decav81B0_M0216_01.nc, global attribute `summary`):
      "Mixed layer is calculated for each profile by estimating the depth for which the
       potential density at 10m (reference depth) increases by 0.125 kg*m-3."

  FESOM MLD2 (src/oce_ale_pressure_bv.F90:462):
      rhopot(nz) - rhopot(nzmin) > sigma_theta_crit,  sigma_theta_crit = 0.125 kg/m3
  where nzmin is the TOP MODEL LEVEL.  Same threshold, same variable (potential density),
  different reference depth: 10 m for WOA, the top level for FESOM -- 2.5 m on this mesh.
  MLD2 is additionally initialised to Z(nzmin+1), so it can never return a value
  shallower than 7.5 m, whereas a 10 m-referenced MLD cannot return less than 10 m.

WHAT THIS SCRIPT DOES.  Recomputes MLD from the model's own monthly T/S under both
conventions, on the native mesh, and reports the difference.  It first reproduces FESOM's
own MLD2 from the same fields as a check that the convention has been read correctly --
if that check fails, nothing below it means anything.

EOS NOTE.  FESOM uses Jackett-McDougall 1995; this uses EOS-80 (the `seawater` package),
which is the convention WOA is built on.  Over the top 100 m the two differ by
O(0.01 kg/m3) against a 0.125 threshold, which is why the reproduction check is expected
to be close but not exact.

MESH GUARD.  core3 = 220509 nodes.  core3_beta (211567) is a different mesh and mixing
them silently produces nonsense.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, seawater as sw, warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESHF = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
ARM = os.environ.get('ARM', '11X')
YEAR = int(os.environ.get('YEAR', '1389'))
CRIT = 0.125          # kg/m3, identical in both conventions
REF_WOA = 10.0        # m, WOA reference depth

BANDS = [('90-60S  Antarctic/SO', -90, -60), ('60-45S  subantarctic', -60, -45),
         ('45-30S', -45, -30), ('30S-30N  tropics', -30, 30), ('30-45N', 30, 45),
         ('45-60N  subpolar NA', 45, 60), ('60-90N  Arctic/Nordic', 60, 90)]

mesh = xr.open_dataset(MESHF)
NOD = mesh.sizes['ncells'] if 'ncells' in mesh.sizes else mesh['lon'].size
assert NOD == 220509, f'expected core3 (220509 nodes), got {NOD} -- wrong mesh'
lat = mesh['lat'].values
area = mesh['cell_area'].values
cav = mesh['cav_nod_mask'].values.astype(bool)

d = f'{R}/{ARM}/outdata/fesom'
T = xr.open_dataset(f'{d}/temp.fesom.{YEAR}.nc', decode_times=False)
S = xr.open_dataset(f'{d}/salt.fesom.{YEAR}.nc', decode_times=False)
M2 = xr.open_dataset(f'{d}/MLD2.fesom.{YEAR}.nc', decode_times=False)
z = T['nz'].values.astype(float)                      # level mid-depths, positive down
print(f'{ARM} {YEAR}: {NOD} nodes, {z.size} levels, top three at {z[:3]} m')
print(f'cavity nodes excluded: {cav.sum()}')


def mld_from_profile(sig, zlev, ref_sigma, first_level):
    """Depth where sig exceeds ref_sigma by CRIT, linearly interpolated.

    sig          (nod, nz) potential density anomaly, NaN below the bottom
    ref_sigma    (nod,)    reference density
    first_level  index of the shallowest level allowed to be returned; mirrors FESOM
                 initialising MLD2 to Z(nzmin+1) rather than to the surface.
    Nodes whose profile never crosses the threshold get the deepest valid depth, which
    is the physically meaningful answer (mixed to the bottom), not a missing value.
    """
    dsig = sig - ref_sigma[:, None]
    valid = np.isfinite(dsig)
    over = valid & (dsig > CRIT)
    over[:, :first_level] = False
    any_over = over.any(axis=1)
    k = np.argmax(over, axis=1)                       # first True; 0 where none
    k = np.maximum(k, first_level)
    km1 = k - 1
    n = np.arange(sig.shape[0])
    d0, d1 = dsig[n, km1], dsig[n, k]
    z0, z1 = zlev[km1], zlev[k]
    with np.errstate(invalid='ignore', divide='ignore'):
        out = z0 + (CRIT - d0) * (z1 - z0) / (d1 - d0)
    out = np.where(np.isfinite(out), out, zlev[k])
    # profiles that never cross: mixed to the bottom
    deepest = np.where(valid.any(axis=1), zlev[valid.shape[1] - 1 -
                       np.argmax(valid[:, ::-1], axis=1)], np.nan)
    out = np.where(any_over, out, deepest)
    return np.where(valid.any(axis=1), out, np.nan)


acc = {k: np.zeros(NOD) for k in ('fesom_conv', 'woa_conv', 'model_mld2')}
cnt = np.zeros(NOD)
for it in range(T.sizes['time']):
    t = T['temp'].values[it].astype(np.float64)       # (nod, nz)
    s = S['salt'].values[it].astype(np.float64)
    bad = ~np.isfinite(t) | ~np.isfinite(s) | (s <= 0)
    t[bad] = np.nan; s[bad] = np.nan
    sig = sw.dens0(s, t) - 1000.0                     # sigma-theta, EOS-80

    # FESOM's convention: reference = top level, floor at level 1
    ref_f = sig[:, 0]
    mld_f = mld_from_profile(sig, z, ref_f, first_level=1)

    # WOA's convention: reference interpolated to 10 m, floor at the first level below 10 m
    k10 = int(np.searchsorted(z, REF_WOA))            # first level deeper than 10 m
    w = (REF_WOA - z[k10 - 1]) / (z[k10] - z[k10 - 1])
    ref_w = sig[:, k10 - 1] + w * (sig[:, k10] - sig[:, k10 - 1])
    ref_w = np.where(np.isfinite(ref_w), ref_w, sig[:, k10 - 1])
    mld_w = mld_from_profile(sig, z, ref_w, first_level=k10)

    m2 = np.abs(M2['MLD2'].values[it].astype(np.float64))
    ok = np.isfinite(mld_f) & np.isfinite(mld_w) & np.isfinite(m2) & ~cav
    acc['fesom_conv'][ok] += mld_f[ok]
    acc['woa_conv'][ok] += mld_w[ok]
    acc['model_mld2'][ok] += m2[ok]
    cnt[ok] += 1

good = cnt > 0
for k in acc:
    acc[k] = np.where(good, acc[k] / np.maximum(cnt, 1), np.nan)

print('\n=== CHECK: does the offline recomputation reproduce FESOM MLD2? ===')
dd = acc['fesom_conv'][good] - acc['model_mld2'][good]
w = area[good]
print(f'  offline(FESOM convention) - model MLD2:  mean {np.average(dd, weights=w):+.2f} m, '
      f'median {np.median(dd):+.2f} m, p5 {np.percentile(dd,5):+.2f}, p95 {np.percentile(dd,95):+.2f}')
print('  (residual is the EOS-80 vs Jackett-McDougall difference plus interpolation detail)')

print('\n=== THE DEFINITIONAL OFFSET: 10 m reference minus top-level reference ===')
print(f'{"band":<24}{"n":>8}{"FESOM conv":>12}{"WOA conv":>11}{"diff":>9}{"diff %":>9}')
for name, la, lb in BANDS:
    sel = good & (lat >= la) & (lat < lb)
    if sel.sum() == 0:
        continue
    w = area[sel]
    a = np.average(acc['fesom_conv'][sel], weights=w)
    b = np.average(acc['woa_conv'][sel], weights=w)
    print(f'{name:<24}{sel.sum():>8}{a:>12.2f}{b:>11.2f}{b-a:>9.2f}{100*(b-a)/a:>8.1f}%')
sel = good
w = area[sel]
a = np.average(acc['fesom_conv'][sel], weights=w)
b = np.average(acc['woa_conv'][sel], weights=w)
print(f'{"GLOBAL":<24}{sel.sum():>8}{a:>12.2f}{b:>11.2f}{b-a:>9.2f}{100*(b-a)/a:>8.1f}%')

print('\n=== how often does the floor bind? ===')
for label, arr, floor in (('FESOM conv (7.5 m)', acc['fesom_conv'], z[1]),
                          ('WOA conv (15 m)', acc['woa_conv'], z[int(np.searchsorted(z, REF_WOA))])):
    at = good & (arr < floor + 1e-6)
    print(f'  {label:<22} {100*area[at].sum()/area[good].sum():5.2f}% of ocean area sits at the floor')
np.savez(os.path.join(os.path.dirname(__file__), f'mld_audit_{ARM}_{YEAR}.npz'),
         lat=lat, area=area, good=good, **acc)
print('\nsaved per-node fields alongside this script')
