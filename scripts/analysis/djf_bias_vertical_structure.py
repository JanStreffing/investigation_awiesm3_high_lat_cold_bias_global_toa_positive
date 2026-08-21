"""Is the DJF land cold bias SURFACE-CONFINED, or is the whole column cold?

WHY THIS DECIDES THE NEXT LEVER.  land_bias_season_band.py localised the coupled land
cold bias to winter (land-minus-ocean gap -1.37 K DJF against -0.09 K JJA).  Every
remaining candidate -- stable-regime mixing (vdfexcu ZCB/ZLMIN), thermal roughness
(ECE_TUNE_RVZ0H), skin conductivity -- works by REDISTRIBUTING heat within the column:
more coupling warms the screen and cools the air just above it.  Such a lever can only
work if the bias is concentrated at the surface.  If the whole lower troposphere is
equally cold, mixing moves the error around instead of removing it, and the cause is
radiative or advective instead.

Two levers were scored null today on exactly this bias -- ECE_LAMSK_SN (skin
conductivity, lamsk_j_series_band.py) and Raupach z0 (raupach_z0_band.py).  Rather than
try a third by hunch, measure whether the redistribution class can work at all.

METRIC.  Model minus ERA5 at 2 m, 1000, 925 and 850 hPa, plus the low-level inversion
strength T925 - T2m and the surface inversion T2m - Tskin, over NH land by band.

  * bias shrinking upward  -> surface-confined; the inversion is too strong; a coupling
    or mixing lever is the right class, and the target is the inversion-strength error.
  * bias roughly uniform   -> deep cold; redistribution is the wrong class entirely.

REFERENCE.  ERA5 monthly T on pressure levels from /pool/data/ERA5/E5/pl/an/1M/130,
DJF 1990-2014, expanded from N320 reduced Gaussian and bilinearly remapped to the model
grid (see ERA5_T_pl_DJF_1990-2014.nc).  1990-2014 matches 11P's constant-1990 forcing.
ERA5 is an HTESSEL sibling, so it is suggestive rather than authoritative for a land
surface claim -- but the vertical STRUCTURE of a 2-3 K bias is far outside the range
where that caveat bites.

Screen level against 925 hPa is compared as a difference of differences, so the
systematic part of any model-ERA5 offset cancels out of the inversion metric.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
E5PL = '/work/ab0246/a270092/obs/era5/netcdf/ERA5_T_DJF_%d.nc'
E5T2M = '/work/ab0246/a270092/obs/era5/netcdf/T2M.nc'
E5SKT = '/work/ab0246/a270092/obs/era5/netcdf/SKT_mon.nc'  # 1979-2008
Y0, Y1 = 1380, 1389
DJF = [11, 0, 1]                      # month indices for Dec, Jan, Feb
ARMS = [('11G inppmin50k', f'{R}/Tuning_test_11G_inppmin50k'),
        ('11N lx4 1850', f'{R}/11N'),
        ('11P lx4 1990', f'{R}/11P')]
BANDS = [('60-90N', 60, 90), ('30-60N', 30, 60), ('ALL LAND', -90, 90)]


def grid():
    f = sorted(glob.glob(f'{R}/11P/outdata/oifs/atm_remapped_1m_lsm_*.nc'))[0]
    with xr.open_dataset(f, decode_times=False) as d:
        m = np.squeeze(d['lsm'].values)
        if m.ndim == 3:
            m = m[0]
        return m, np.squeeze(d['lat'].values)


msk, lat = grid()
land = msk > 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], msk.shape).copy()


def am(f, sel):
    m = sel & np.isfinite(f)
    return float(np.average(f[m], weights=W[m])) if m.any() else np.nan


def model_djf(path, var, plev=None):
    """DJF-mean field over Y0..Y1, optionally at a pressure level [Pa]."""
    acc = []
    for y in range(Y0, Y1 + 1):
        f = f'{path}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(f):
            return None
        with xr.open_dataset(f, decode_times=False) as d:
            k = 't' if plev is not None else [c for c in d.data_vars
                                              if 'bounds' not in c][0]
            a = d[k]
            if plev is not None:
                pl = np.squeeze(d['pressure_levels'].values)
                a = a.isel({[x for x in a.dims if 'press' in x][0]:
                            int(np.argmin(np.abs(pl - plev)))})
            a = np.squeeze(a.values)
        acc.append(a[DJF].mean(0))
    a = np.mean(acc, axis=0)
    return a - 273.15 if np.nanmean(a) > 100 else a


def era5_pl(plev):
    """One file per level.  A merged multi-level file is NOT used: cdo wrote it with
    dims (time, lat, plev, lon) and slicing it by name returned a near-constant field
    (251-253 K globally at 850 hPa), so single-level files remove the ambiguity."""
    with xr.open_dataset(E5PL % plev, decode_times=False) as d:
        k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        a = np.squeeze(d[k].values)
    if a.shape != msk.shape:
        raise SystemExit(f'ERA5 {plev} Pa has shape {a.shape}, expected {msk.shape}')
    return a - 273.15 if np.nanmean(a) > 100 else a


def below_ground(path, plev):
    """True where the pressure level lies UNDER the surface, so both model and ERA5
    are reporting downward extrapolation rather than air.  Over land these cells are
    the difference between a +28 K 'bias' and a real one."""
    acc = []
    for y in range(Y0, Y1 + 1):
        f = f'{path}/outdata/oifs/atm_remapped_1m_sp_{y}-{y}.nc'
        with xr.open_dataset(f, decode_times=False) as d:
            a = np.squeeze(d['sp'].values)
        acc.append(a[DJF].mean(0))
    return np.mean(acc, axis=0) < plev


def era5_t2m():
    """DJF climatology.  Assert the month axis before trusting it -- averaging all
    300 steps into one annual field once produced a fictitious -22 K DJF bias."""
    with xr.open_dataset(E5T2M, decode_times=True) as d:
        clim = d['t2m'].groupby('time.month').mean('time')
        la = [c for c in clim.dims if 'lat' in c][0]
        lo = [c for c in clim.dims if 'lon' in c][0]
        out = np.asarray(clim.interp({la: ('lat', lat),
                                      lo: ('lon', np.linspace(0.45, 359.55, msk.shape[1]))}
                                     ).values, dtype=float)
    out = out - 273.15 if np.nanmean(out) > 100 else out
    nh = lat > 30
    cyc = [float(np.average(out[i][nh].mean(axis=1), weights=np.cos(np.deg2rad(lat))[nh]))
           for i in range(12)]
    if int(np.argmax(cyc)) not in (5, 6, 7):
        raise SystemExit(f'ERA5 month axis wrong: NH peak at index {np.argmax(cyc)}')
    return out[DJF].mean(0)


def era5_skt():
    """DJF skin-temperature climatology on the model grid.  Only ever used inside the
    difference-of-differences (T2m - Tskt), so the 1979-2008 window against the pressure
    levels' 1990-2014 does not bite: a period offset moves both terms together."""
    with xr.open_dataset(E5SKT, decode_times=True) as d:
        clim = d['skt'].groupby('time.month').mean('time')
        la = [c for c in clim.dims if 'lat' in c][0]
        lo = [c for c in clim.dims if 'lon' in c][0]
        out = np.asarray(clim.interp(
            {la: ('lat', lat),
             lo: ('lon', np.linspace(0.45, 359.55, msk.shape[1]))}).values, dtype=float)
    out = out - 273.15 if np.nanmean(out) > 100 else out
    return out[DJF].mean(0)


o2, o1000, o925, o850 = era5_t2m(), era5_pl(100000), era5_pl(92500), era5_pl(85000)
oskt = era5_skt()

print('DJF model-minus-ERA5 bias by height, NH land.  ERA5 pl 1990-2014; '
      f'model {Y0}-{Y1}.\n')
hdr = ['T2m', 'T1000', 'T925', 'T850', 'lowinv bias', 'sfcinv bias']
for label, path in ARMS:
    m2 = model_djf(path, '2t')
    if m2 is None:
        print(f'{label}: incomplete\n')
        continue
    mskt = model_djf(path, 'skt')
    m1000 = model_djf(path, 'pl_t', 100000)
    m925 = model_djf(path, 'pl_t', 92500)
    m850 = model_djf(path, 'pl_t', 85000)
    bg = {p: below_ground(path, p) for p in (100000, 92500, 85000)}
    print(f'=== {label} ===   below-ground land cells excluded per level: '
          + ', '.join(f'{p//100} hPa {int((bg[p] & land).sum())}'
                      for p in (100000, 92500, 85000)))
    print(f'{"band":>9} | ' + ' | '.join(f'{h:>16}' for h in hdr))
    for bname, lo, hi in BANDS:
        sel = land & np.broadcast_to(((lat >= lo) & (lat < hi))[:, None], msk.shape)
        b2 = am(m2 - o2, sel)
        s1000, s925, s850 = (sel & ~bg[100000], sel & ~bg[92500], sel & ~bg[85000])
        vals = [b2, am(m1000 - o1000, s1000), am(m925 - o925, s925),
                am(m850 - o850, s850),
                am((m925 - m2) - (o925 - o2), s925),
                am((m2 - mskt) - (o2 - oskt), sel)]
        print(f'{bname:>9} | ' + ' | '.join(f'{v:+16.3f}' for v in vals))
    print()
print('Both inversion columns are BIASES (model minus ERA5), positive = model '
      'inversion too strong.\nlowinv = T925-T2m (air-side, aerodynamic); '
      'sfcinv = T2m-Tskt (skin decoupling).')
