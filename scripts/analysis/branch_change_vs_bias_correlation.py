"""Does a branch's T2m change undo the control's T2m bias?  Pattern correlation and projection.

Bias = CONTROL minus ERA5 1990-2014, change = branch minus CONTROL, both means over Y0-Y1 on
the model's remapped grid.  A lever that corrects the bias has a change pattern
anti-correlated with the bias.  Reported per season and domain, area-weighted:
  r_grid  centred pattern correlation on grid cells
  r_box   the same on 10x10 degree box means (suppresses grid-scale noise)
  slope   regression of change on centred bias: -1 would remove the pattern exactly,
          0 leaves it untouched
  frac    fraction of the bias variance removed, 1 - var(bias+change)/var(bias),
          centred, so a uniform shift does not count
Centring removes the domain mean from both fields: a pre-industrial run should be uniformly
colder than present-day ERA5, so the useful question is about the pattern, and the domain
means are listed separately.

Usage:  CONTROL=PI200 RUNS=PI200_spp,PI200_mle,PI200_h0,PI200_all3 Y0=1565 Y1=1579 \
        python3 scripts/analysis/branch_change_vs_bias_correlation.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
CONTROL = os.environ.get('CONTROL', 'PI200')
RUNS = os.environ.get('RUNS', 'PI200_spp,PI200_mle,PI200_h0,PI200_all3').split(',')
Y0, Y1 = int(os.environ.get('Y0', 1565)), int(os.environ.get('Y1', 1579))
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'
SEAS = {'ANN': list(range(12)), 'DJF': [11, 0, 1], 'JJA': [5, 6, 7]}
DOMAINS = [('global', -90, 90, 'all'), ('45-90S', -90, -45, 'all'), ('45-90S ocean', -90, -45, 'ocean'),
           ('60-90S', -90, -60, 'all'), ('45-90N', 45, 90, 'all')]


def load(arm):
    acc = 0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{arm}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc', decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            acc = acc + np.squeeze(d[k].values).astype('f8') / (Y1 - Y0 + 1)
            lat = np.squeeze(d['lat'].values); lon = np.squeeze(d['lon'].values)
    return acc, lat, lon


C, lat, lon = load(CONTROL)
order = np.argsort(lon % 360); lon = (lon % 360)[order]; C = C[..., order]
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = (lsm[0] if lsm.ndim == 3 else lsm)[:, order]
with xr.open_dataset('/work/ab0246/a270092/obs/era5/netcdf/T2M.nc') as d:
    c = d['t2m'].groupby('time.month').mean('time')
    la = [x for x in c.dims if 'lat' in x][0]; lo = [x for x in c.dims if 'lon' in x][0]
    c = c.assign_coords({lo: c[lo] % 360}).sortby(lo)
    E = np.asarray(c.interp({la: ('lat', lat), lo: ('lon', lon)}).values, float)
LAT = np.broadcast_to(lat[:, None], lsm.shape); W = np.cos(np.deg2rad(LAT))
iby = np.floor((LAT + 90) / 10).astype(int); ibx = np.floor(np.broadcast_to(lon[None, :], lsm.shape) / 10).astype(int)


def stats(b, d, k):
    k = k & np.isfinite(b) & np.isfinite(d)
    w = W[k]; bb = b[k] - np.average(b[k], weights=w); dd = d[k] - np.average(d[k], weights=w)
    r = np.sum(w * bb * dd) / np.sqrt(np.sum(w * bb ** 2) * np.sum(w * dd ** 2))
    slope = np.sum(w * bb * dd) / np.sum(w * bb ** 2)
    frac = 1 - np.sum(w * (bb + dd) ** 2) / np.sum(w * bb ** 2)
    ids = iby[k] * 36 + ibx[k]                                  # 10-degree boxes
    sw = np.bincount(ids, w); m = sw > 0
    bx = np.bincount(ids, w * b[k])[m] / sw[m]; dx = np.bincount(ids, w * d[k])[m] / sw[m]; wx = sw[m]
    bx -= np.average(bx, weights=wx); dx -= np.average(dx, weights=wx)
    rb = np.sum(wx * bx * dx) / np.sqrt(np.sum(wx * bx ** 2) * np.sum(wx * dx ** 2))
    return r, rb, slope, frac, np.average(b[k], weights=w), np.average(d[k], weights=w)


B = {a: load(a)[0][..., order] for a in RUNS}
print(f'{Y0}-{Y1}: change (branch - {CONTROL}) against bias ({CONTROL} - ERA5).  Negative r = corrects the pattern.')
for s, ms in SEAS.items():
    bias = C[ms].mean(0) - E[ms].mean(0)
    print(f'\n{s}')
    print(f'  {"domain":<14}{"run":<12}{"r_grid":>8}{"r_box":>8}{"slope":>8}{"frac":>8}{"mean bias":>11}{"mean chg":>10}')
    for name, a, b, surf in DOMAINS:
        k = (LAT >= a) & (LAT <= b)
        if surf == 'ocean':
            k &= lsm <= 0.5
        for arm in RUNS:
            chg = B[arm][ms].mean(0) - C[ms].mean(0)
            r, rb, sl, fr, mb, mc = stats(bias, chg, k)
            print(f'  {name:<14}{arm.replace(CONTROL + "_", ""):<12}{r:+8.2f}{rb:+8.2f}{sl:+8.2f}{fr:+8.2f}{mb:+11.2f}{mc:+10.2f}')
