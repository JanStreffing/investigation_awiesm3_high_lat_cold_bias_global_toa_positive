"""Base -> last best -> the SB lever, all at 1850 forcing, in one table.

THE THREE ARMS, all PI/1850 so the comparison between them is clean:
  11E  Tuning_test_11E_swemin15_K1  the campaign base
  11N  LX4 stack coupled            the best run before this round
  11Q  11N + RSBLB 2.0              stable-BL mixing, this round

READ THE COLUMNS, NOT THE ABSOLUTE NUMBERS, FOR ANYTHING vs ERA5.  All three run 1850
forcing and ERA5 is present-day, so every T2m bias here carries a forcing offset of order
a kelvin that is NOT model error.  Only 11P/11R are period-matched.  The 11E->11N->11Q
DIFFERENCES are the meaningful quantity, and they are what this table is for.  CERES and
sea-ice comparisons are less affected but not immune.

Common window 1360-1389 (30 yr), the first decade discarded as spin-up.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
ARMS = [('11E base', f'{R}/Tuning_test_11E_swemin15_K1'),
        ('11N LX4', f'{R}/11N'),
        ('11Q +RSBLB2', f'{R}/11Q')]
E = '/work/ab0246/a270092/obs/era5/netcdf'
Y0, Y1 = 1360, 1389
DJF, JJA = [11, 0, 1], [5, 6, 7]

f = sorted(glob.glob(f'{R}/11N/outdata/oifs/atm_remapped_1m_lsm_*.nc'))[0]
with xr.open_dataset(f, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values)
    m = m[0] if m.ndim == 3 else m
    lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
land, ocean = m > 0.5, m <= 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], m.shape).copy()
bnd = lambda lo, hi: np.broadcast_to(((lat >= lo) & (lat < hi))[:, None], m.shape)
arctic, mid = land & bnd(60, 90), land & bnd(30, 60)
sib = land & bnd(50, 70) & np.broadcast_to(((lon >= 60) & (lon <= 140))[None, :], m.shape)
trop, so = bnd(-30, 30), bnd(-65, -45)
glob_ = np.ones_like(m, dtype=bool)
RE = 6.371e6
AREA = np.broadcast_to((RE**2 * np.cos(np.deg2rad(lat)) * np.deg2rad(abs(lat[1]-lat[0]))
                        * 2*np.pi/m.shape[1])[:, None], m.shape)


def am(f_, s):
    k = s & np.isfinite(f_)
    return float(np.average(f_[k], weights=W[k])) if k.any() else np.nan


def load(path, var, plev=None):
    out = []
    for y in range(Y0, Y1 + 1):
        p = f'{path}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(p):
            return None
        with xr.open_dataset(p, decode_times=False) as d:
            k = 't' if plev is not None else \
                [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = d[k]
            if plev is not None:
                pl = np.squeeze(d['pressure_levels'].values)
                a = a.isel({[x for x in a.dims if 'press' in x][0]:
                            int(np.argmin(np.abs(pl - plev)))})
            out.append(np.squeeze(a.values))
    return np.stack(out)


def era_pl(p):
    with xr.open_dataset(f'{E}/ERA5_T_DJF_{p}.nc', decode_times=False) as d:
        a = np.squeeze(d[[k for k in d.data_vars if 'bnds' not in k][0]].values)
    return a - 273.15 if np.nanmean(a) > 100 else a


with xr.open_dataset(f'{E}/T2M.nc') as d:
    c = d['t2m'].groupby('time.month').mean('time')
    la = [x for x in c.dims if 'lat' in x][0]
    lo = [x for x in c.dims if 'lon' in x][0]
    O = np.asarray(c.interp({la: ('lat', lat),
                             lo: ('lon', np.linspace(0.45, 359.55, m.shape[1]))}).values, float)
O = O - 273.15 if np.nanmean(O) > 100 else O
nh = lat > 30
cyc = [float(np.average(O[i][nh].mean(1), weights=np.cos(np.deg2rad(lat))[nh])) for i in range(12)]
assert int(np.argmax(cyc)) in (5, 6, 7), 'ERA5 month axis wrong'
o2_djf, o2_jja, o2_ann = O[DJF].mean(0), O[JJA].mean(0), O.mean(0)
o925 = era_pl(92500)

rows = {}
for label, path in ARMS:
    t2 = load(path, '2t')
    if t2 is None:
        print(f'{label}: incomplete'); continue
    t2 = t2 - 273.15 if np.nanmean(t2) > 100 else t2
    p925 = load(path, 'pl_t', 92500)
    p925 = p925 - 273.15 if np.nanmean(p925) > 100 else p925
    sp = np.mean([load(path, 'sp')[i][DJF].mean(0) for i in range(1)], axis=0) \
        if False else np.mean(load(path, 'sp')[:, DJF].mean(1), axis=0)
    tsr, tsrc, ttr = load(path, 'tsr'), load(path, 'tsrc'), load(path, 'ttr')
    ci = load(path, 'ci')
    n = t2.shape[0]
    d = {}
    md = t2[:, DJF].mean(1).mean(0); mj = t2[:, JJA].mean(1).mean(0); ma_ = t2.mean(1).mean(0)
    for rn, sel in [('60-90N', arctic), ('30-60N', mid), ('all land', land)]:
        d[f'DJF T2m bias {rn}'] = am(md - o2_djf, sel)
    d['ANN T2m bias all land'] = am(ma_ - o2_ann, land)
    d['ANN T2m bias ocean'] = am(ma_ - o2_ann, ocean)
    d['JJA T2m bias Siberia'] = am(mj - o2_jja, sib)
    s9 = arctic & (sp >= 92500)
    d['DJF inversion bias 60-90N'] = am((p925[:, DJF].mean(1).mean(0) - md) - (o925 - o2_djf), s9)
    net = (tsr + ttr).mean(1).mean(0) / 3600.
    d['net TOA global'] = am(net, glob_)
    swcre = (tsr - tsrc).mean(1).mean(0) / 3600.
    d['SW CRE SO 45-65S'] = am(swcre, so)
    d['SW CRE tropics'] = am(swcre, trop)
    d['NH ice Mar [1e6km2]'] = float((ci[:, [2]].mean(1).mean(0) * AREA)[ocean & bnd(45, 90)].sum() / 1e12)
    d['SH ice Sep [1e6km2]'] = float((ci[:, [8]].mean(1).mean(0) * AREA)[ocean & bnd(-90, -45)].sum() / 1e12)
    rows[label] = d

keys = list(next(iter(rows.values())).keys())
labs = list(rows.keys())
print(f'1850-forcing arms, {Y0}-{Y1}.  T2m/inversion vs ERA5; the OFFSET is forcing, '
      f'the DIFFERENCES are the result.\n')
print(f'{"metric":>28} | ' + ' | '.join(f'{l:>12}' for l in labs)
      + ' | ' + f'{"11Q-11N":>9} {"11Q-11E":>9}')
print('-' * 100)
for k in keys:
    v = [rows[l][k] for l in labs]
    print(f'{k:>28} | ' + ' | '.join(f'{x:12.3f}' for x in v)
          + ' | ' + f'{v[2]-v[1]:+9.3f} {v[2]-v[0]:+9.3f}')
