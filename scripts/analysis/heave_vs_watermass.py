"""Is the deep warming carried by moving isopycnals, or by warming the water on them?

THE SPLIT.  ohc_drift_by_depth_band.py shows the model losing heat from 60-340 m and
gaining it at 340-2000 m over 30 years, peaking at +0.50 K near 600 m.  Two mechanisms
produce that at fixed depth and they call for opposite fixes:

  HEAVE       the density surfaces move (deeper thermocline, changed gyre spin-up), and a
              fixed depth simply samples different water.  Temperature ON a density
              surface is unchanged.  This is circulation, not mixing.
  WATER-MASS  the water on a given density surface actually gets warmer -- heat crossed
              density surfaces, or the water was formed with different properties.  This
              is diapycnal mixing or a change in the source water.

The discriminator is theta on isopycnals.  If d(theta)|sigma is near zero while
d(theta)|z is large, it is heave.  If they track each other, heat is crossing isopycnals.

METHOD.  Decade means of theta and S; sigma_0 from EOS-80; theta and depth interpolated
onto a fixed sigma_0 grid per node (sigma made monotonic by running maximum, which is what
a sorted/stably-stratified column would give); band means taken in the mesh's own
geometry, volume-weighted by the layer thickness each sigma level occupies.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, seawater as sw, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESHF = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
ARM = os.environ.get('ARM', '11X')
EARLY = list(range(1350, 1360))
LATE = list(range(1380, 1390))
SIG = np.arange(24.0, 28.05, 0.05)

BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45), ('45-30S', -45, -30),
         ('30S-30N', -30, 30), ('30-45N', 30, 45), ('45-60N', 45, 60), ('60-90N', 60, 90)]


def decade(years, var):
    acc, n = None, 0
    for y in years:
        with xr.open_dataset(f'{R092}/{ARM}/outdata/fesom/{var}.fesom.{y}.nc',
                             decode_times=False) as d:
            a = d[var].values.astype('f8').mean(axis=0)
        acc = a if acc is None else acc + a
        n += 1
    return acc / n


def on_sigma(T, S, zmid):
    """theta(sigma) and z(sigma) per node, NaN where the sigma level is not present."""
    sig = sw.dens0(S, T) - 1000.0
    nnod, nz = T.shape
    th_out = np.full((nnod, SIG.size), np.nan)
    z_out = np.full((nnod, SIG.size), np.nan)
    smon = np.fmax.accumulate(np.where(np.isfinite(sig), sig, -np.inf), axis=1)
    smon[~np.isfinite(sig)] = np.nan
    for i in range(nnod):
        s = smon[i]
        good = np.isfinite(s)
        if good.sum() < 3:
            continue
        sv = s[good]
        tv = T[i][good]
        zv = zmid[good]
        d = np.diff(sv) > 1e-6
        keep = np.concatenate([[True], d])
        sv, tv, zv = sv[keep], tv[keep], zv[keep]
        if sv.size < 3:
            continue
        inside = (SIG >= sv[0]) & (SIG <= sv[-1])
        if not inside.any():
            continue
        th_out[i, inside] = np.interp(SIG[inside], sv, tv)
        z_out[i, inside] = np.interp(SIG[inside], sv, zv)
    return th_out, z_out


def main():
    with xr.open_dataset(MESHF, decode_times=False) as m:
        lat = m['lat'].values.astype('f8')
        area = m['cell_area'].values.astype('f8')
    assert lat.size == 220509, 'MESH GUARD: not core3'
    with xr.open_dataset(f'{R092}/{ARM}/outdata/fesom/temp.fesom.{LATE[-1]}.nc',
                         decode_times=False) as d:
        zmid = d['nz'].values.astype('f8')

    print('loading decade means ...', flush=True)
    Ta, Sa = decade(EARLY, 'temp'), decade(EARLY, 'salt')
    Tb, Sb = decade(LATE, 'temp'), decade(LATE, 'salt')
    print('projecting onto sigma_0 ...', flush=True)
    tha, za = on_sigma(Ta, Sa, zmid)
    thb, zb = on_sigma(Tb, Sb, zmid)

    dth = thb - tha           # water-mass / diapycnal component
    dz = zb - za              # heave: how far the surface moved (+ = deeper)
    zbar = 0.5 * (za + zb)

    print('\n' + '=' * 100)
    print(f'{ARM}: change ON sigma_0 surfaces, {EARLY[0]}-{EARLY[-1]} to {LATE[0]}-{LATE[-1]}')
    print('  d(theta)|sigma  = heat that CROSSED density surfaces (diapycnal / water-mass)')
    print('  d(z)|sigma      = how far the surface moved (heave; + is deeper)')
    print('=' * 100)
    for name, lo, hi in BANDS:
        sel = (lat >= lo) & (lat < hi)
        print(f'\n{name}')
        print(f'{"sigma_0":>9}{"mean depth":>12}{"d(theta)|sig":>15}{"d(z)|sig":>12}')
        for j in range(0, SIG.size, 4):
            g = sel & np.isfinite(dth[:, j]) & np.isfinite(dz[:, j])
            if g.sum() < 200:
                continue
            w = area[g]
            zz = np.sum(zbar[g, j] * w) / w.sum()
            if zz < 80 or zz > 2500:
                continue
            t = np.sum(dth[g, j] * w) / w.sum()
            h = np.sum(dz[g, j] * w) / w.sum()
            print(f'{SIG[j]:>9.2f}{zz:>12.0f}{t:>+15.3f}{h:>+12.1f}')
    np.savez(os.path.join('data', f'heave_watermass_{ARM}.npz'),
             SIG=SIG, dth=dth, dz=dz, zbar=zbar, lat=lat, area=area)
    print('\nsaved data/heave_watermass_11X.npz')


if __name__ == '__main__':
    main()
