"""Nino3.4 mean SST and ENSO amplitude for the coupled arms against HadISST.

Nino3.4 box: 5S-5N, 170W-120W.  Model: FESOM daily `sst` on mesh nodes, box nodes
area-weighted with mesh.nc cell_area, averaged to calendar months.  Per run, over its last
NYR complete years:
  mean      box-mean SST [degC]
  cycle     amplitude of the mean seasonal cycle (max - min of the monthly climatology) [K]
  enso_sd   std of monthly anomalies after removing the monthly climatology and a linear
            trend [K]  (the ENSO amplitude)
HadISST over matched periods: 1980-2009 for the 1990-forced arms, 1870-1899 as the
pre-industrial proxy for the 1850 arms (sparse; a proxy, not a PI truth).

Usage:  python3 scripts/analysis/nino34_enso_state.py        (heavy: run on Slurm)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
from multiprocessing import Pool
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESH = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
HAD = '/work/ab0246/a270092/obs/hadisst2/HadISST_sst.nc'
NYR = int(os.environ.get('NYR', 30))
# run: (forcing, last year)
RUNS = {'11X': ('1850', 1389), '15F': ('1850', 1379), '16C': ('1850', 1389), '16D': ('1850', 1389),
        '16E': ('1850', 1389), '11R': ('1990', 1389), '11V': ('1990', 1389), '16E_1990': ('1990', 1399)}

with xr.open_dataset(MESH) as m:
    AREA = m['cell_area'].values.astype('f8')


def box_nodes(path):
    with xr.open_dataset(path, decode_times=False) as d:
        lat, lon = d['lat'].values, d['lon'].values % 360
    k = (lat >= -5) & (lat <= 5) & (lon >= 190) & (lon <= 240)
    return np.where(k)[0]


def year_months(args):
    run, y, idx = args
    p = f'{R}/{run}/outdata/fesom/sst.fesom.{y}.nc'
    if not os.path.exists(p):
        return None
    with xr.open_dataset(p) as d:
        s = d['sst'].isel(nod2=idx).values.astype('f8')          # (days, nodes)
        month = np.array([int(str(x)[5:7]) for x in d['time'].values])   # model years predate pandas' datetime range
    if month.size < 365:
        return None
    w = AREA[idx] if AREA.size >= idx.max() + 1 else np.ones(idx.size)
    daily = (s * w).sum(1) / w.sum()
    return y, [float(daily[month == mo].mean()) for mo in range(1, 13)]


def stats(monthly):
    a = np.asarray(monthly)                                         # (years, 12)
    clim = a.mean(0); anom = (a - clim).ravel()
    x = np.arange(anom.size); anom = anom - np.polyval(np.polyfit(x, anom, 1), x)
    return a.mean(), clim.max() - clim.min(), anom.std()


def hadisst(y0, y1):
    with xr.open_dataset(HAD) as d:
        s = d['sst'].sel(latitude=slice(5, -5), longitude=slice(-170, -120))
        s = s.sel(time=slice(f'{y0}-01-01', f'{y1}-12-31'))
        s = s.where(s > -100)
        w = np.cos(np.deg2rad(s['latitude']))
        ts = s.weighted(w).mean(('latitude', 'longitude')).values
    return stats(ts[: (ts.size // 12) * 12].reshape(-1, 12))


if __name__ == '__main__':
    jobs = []
    for run, (forc, y1) in RUNS.items():
        f0 = f'{R}/{run}/outdata/fesom/sst.fesom.{y1}.nc'
        if not os.path.exists(f0):
            print(f'{run}: no sst for {y1}, skipped'); continue
        idx = box_nodes(f0)
        jobs += [(run, y, idx) for y in range(y1 - NYR + 1, y1 + 1)]
    with Pool(int(os.environ.get('NPROC', 32))) as pool:
        res = pool.map(year_months, [(r, y, i) for r, y, i in jobs])
    by = {}
    for (run, y, _), r in zip(jobs, res):
        if r:
            by.setdefault(run, {})[r[0]] = r[1]
    print(f'Nino3.4 (5S-5N, 170W-120W), last {NYR} complete years per run')
    print(f'{"run":10s} {"forcing":8s} {"years":10s} {"mean SST":>9s} {"cycle":>6s} {"ENSO sd":>8s}')
    for run, (forc, y1) in RUNS.items():
        if run not in by:
            continue
        ys = sorted(by[run]); m, c, sd = stats([by[run][y] for y in ys])
        print(f'{run:10s} {forc:8s} {ys[0]}-{ys[-1]} {m:9.2f} {c:6.2f} {sd:8.3f}')
    for lab, (y0, y1) in (('HadISST 1980-2009', (1980, 2009)), ('HadISST 1870-1899', (1870, 1899))):
        m, c, sd = hadisst(y0, y1)
        print(f'{lab:28s}  {m:9.2f} {c:6.2f} {sd:8.3f}')
