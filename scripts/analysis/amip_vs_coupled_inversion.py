"""Does the AMIP base share the coupled winter inversion defect, and by how much?

WHY.  The +2.81 K low-level inversion bias that motivates the SB series (RSBLB) was
measured on COUPLED 11P.  SB0/SB1/SB2 are AMIP arms, so the screen is only meaningful if
the AMIP base carries the same defect.  Checked before spending the compute rather than
after -- and the answer also bounds what an AMIP result can ever prove.

RESULT (2026-08-21), DJF, model minus ERA5, NH land:

    band      | AMIP LX4 lowinv | coupled 11P lowinv
    60-90N    |     +1.54       |      +2.81
    30-60N    |     +0.37       |      +1.07
    ALL LAND  |     +0.43       |      +0.86

Same sign and same vertical shape, so the screen is valid.  But AMIP carries only ~55 %
of the coupled inversion error at 60-90N, and the split is diagnostic: T925 is IDENTICAL
in the two (-0.88), so the whole difference sits at the screen.  ~1.3 K of the coupled
excess therefore has no AMIP counterpart, and the obvious home for it is sea ice, which
AMIP prescribes.  A clean AMIP result cannot close the coupled gap by itself.

CAVEATS.  The AMIP arms are PI-forced and compared against a present-day ERA5, so the
ABSOLUTE T2m offset is expected and is not the quantity read here.  The inversion bias is
a difference of differences, (T925-T2m)_model - (T925-T2m)_ERA5, so a uniform offset
cancels out of it -- that is the whole reason the metric is built that way.

The below-ground mask is taken from 11P's DJF surface pressure: LX4 does not output sp,
and the mask is a topography question, not a simulation one, on an identical grid.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

AMIP = '/work/bb1469/a270092/runtime/oifsamip-cy48/LX4/outdata/oifs'
CPL = '/work/bb1469/a270092/runtime/awiesm3-v3.4/11P/outdata/oifs'
E = '/work/ab0246/a270092/obs/era5/netcdf'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_'
        'ENTSTPC3_CRUNCEPinit/outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
DJF = [11, 0, 1]
BANDS = [('60-90N', 60, 90), ('30-60N', 30, 60), ('ALL LAND', -90, 90)]

with xr.open_dataset(LSMF, decode_times=False) as d:
    msk = np.squeeze(d['lsm'].values)
    msk = msk[0] if msk.ndim == 3 else msk
    lat = np.squeeze(d['lat'].values)
land = msk > 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], msk.shape).copy()


def am(f, s):
    k = s & np.isfinite(f)
    return float(np.average(f[k], weights=W[k]))


def model(path, years, var, plev=None, tag=''):
    acc = []
    for y in years:
        p = f'{path}/atm_remapped_1m_{var}_1m{tag}_{y}-{y}.nc' if tag or path == AMIP \
            else f'{path}/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(p):
            return None
        with xr.open_dataset(p, decode_times=False) as d:
            key = 't' if plev is not None else \
                [k for k in d.data_vars if 'bnds' not in k and 'bounds' not in k][0]
            a = d[key]
            if plev is not None:
                pl = np.squeeze(d['pressure_levels'].values)
                dim = [x for x in a.dims if 'press' in x][0]
                a = a.isel({dim: int(np.argmin(np.abs(pl - plev)))})
            a = np.squeeze(a.values)
        acc.append(a[DJF].mean(0))
    a = np.mean(acc, axis=0)
    return a - 273.15 if np.nanmean(a) > 100 else a


def era_pl(p):
    with xr.open_dataset(f'{E}/ERA5_T_DJF_{p}.nc', decode_times=False) as d:
        a = np.squeeze(d[[k for k in d.data_vars if 'bnds' not in k][0]].values)
    return a - 273.15 if np.nanmean(a) > 100 else a


def era_t2m():
    with xr.open_dataset(f'{E}/T2M.nc') as d:
        c = d['t2m'].groupby('time.month').mean('time')
        la = [x for x in c.dims if 'lat' in x][0]
        lo = [x for x in c.dims if 'lon' in x][0]
        o = np.asarray(c.interp({la: ('lat', lat),
                                 lo: ('lon', np.linspace(0.45, 359.55, msk.shape[1]))
                                 }).values, dtype=float)
    o = o - 273.15 if np.nanmean(o) > 100 else o
    nh = lat > 30
    cyc = [float(np.average(o[i][nh].mean(1), weights=np.cos(np.deg2rad(lat))[nh]))
           for i in range(12)]
    if int(np.argmax(cyc)) not in (5, 6, 7):
        raise SystemExit(f'ERA5 month axis wrong: peak at {np.argmax(cyc)}')
    return o[DJF].mean(0)


sp = np.mean([np.squeeze(xr.open_dataset(f'{CPL}/atm_remapped_1m_sp_{y}-{y}.nc',
                                         decode_times=False)['sp'].values)[DJF].mean(0)
              for y in range(1380, 1390)], axis=0)

o2, o925, o850 = era_t2m(), era_pl(92500), era_pl(85000)
amip_yrs = [y for y in sorted({int(f.split('_')[-1].split('-')[0])
                               for f in glob.glob(f'{AMIP}/atm_remapped_1m_2t_1m_*.nc')})
            if 1872 <= y <= 1915]

for label, path, years, tag in [('AMIP LX4', AMIP, amip_yrs, '_pl'),
                                ('coupled 11P', CPL, list(range(1380, 1390)), '')]:
    m2 = model(path, years, '2t')
    m925 = model(path, years, 'pl_t', 92500, tag)
    m850 = model(path, years, 'pl_t', 85000, tag)
    if m2 is None or m925 is None:
        print(f'{label}: incomplete\n')
        continue
    print(f'=== {label}, DJF {years[0]}-{years[-1]} ({len(years)} yr) ===')
    print(f'{"band":>9} | {"T2m":>7} {"T925":>7} {"T850":>7} | {"lowinv bias":>11}')
    for nm, lo, hi in BANDS:
        s = land & np.broadcast_to(((lat >= lo) & (lat < hi))[:, None], msk.shape)
        s9, s8 = s & (sp >= 92500), s & (sp >= 85000)
        print(f'{nm:>9} | {am(m2 - o2, s):+7.2f} {am(m925 - o925, s9):+7.2f} '
              f'{am(m850 - o850, s8):+7.2f} | {am((m925 - m2) - (o925 - o2), s9):+11.2f}')
    print()
