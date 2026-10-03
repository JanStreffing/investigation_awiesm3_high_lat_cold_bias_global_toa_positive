"""Is the equatorial Pacific cold tongue too cloudy?  Cloud radiative effect against CERES.

Compares what the clouds DO (cloud radiative effect, CRE = all-sky minus clear-sky flux),
which is defined identically in model and satellite, rather than cloud amount: CERES
cldarea is a MODIS mask that counts optically thin cloud and is not comparable to the
model's tcc (printed for reference only).
  TOA  SW CRE = tsr - tsrc,  LW CRE = ttr - ttrc
  SFC  SW CRE = ssr - ssrc,  LW CRE = str - strc   (net, positive down, what the ocean gets)
  SFC  SW down all-sky = ssrd
IFS fluxes are accumulated J/m2 per hourly output step: divided by 3600.
Regions: Nino3.4 (5S-5N, 170W-120W) and the east-Pacific cold tongue (5S-5N, 150W-90W).
CERES EBAF Ed4.1 is a 2000s climatology, so only the 1990-forced arms are compared with
it; the 1850 arm is shown for context.

Usage:  python3 scripts/analysis/nino34_cloud_cre.py
        RUNS=16E_1990:1390:1399,11V:1380:1389 python3 ...
"""
import os
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CER = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
ACC = 3600.0
RUNS = [r.split(':') for r in os.environ.get('RUNS', '16E_1990:1390:1399,11V:1380:1389,16E:1380:1389').split(',')]
BOXES = {'Nino3.4': (-5, 5, 190, 240), 'cold tongue': (-5, 5, 210, 270)}


def mean_field(run, var, y0, y1):
    acc = None; n = 0
    for y in range(int(y0), int(y1) + 1):
        p = f'{R}/{run}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        with xr.open_dataset(p, decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8').mean(0)
            lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values) % 360
        acc = a if acc is None else acc + a; n += 1
    return acc / n, lat, lon


def box_mean(f, lat, lon, box):
    la0, la1, lo0, lo1 = box
    L = np.broadcast_to(lat[:, None], f.shape); O = np.broadcast_to(lon[None, :], f.shape)
    k = (L >= la0) & (L <= la1) & (O >= lo0) & (O <= lo1) & np.isfinite(f)
    return float(np.average(f[k], weights=np.cos(np.deg2rad(L[k]))))


rows = {}
for run, y0, y1 in RUNS:
    F = {v: mean_field(run, v, y0, y1) for v in ('tsr', 'tsrc', 'ttr', 'ttrc', 'ssr', 'ssrc', 'str', 'strc', 'ssrd', 'tcc', 'lcc')}
    lat, lon = F['tsr'][1], F['tsr'][2]
    g = lambda v: F[v][0]
    fields = {'TOA SW CRE': (g('tsr') - g('tsrc')) / ACC, 'TOA LW CRE': (g('ttr') - g('ttrc')) / ACC,
              'SFC SW CRE': (g('ssr') - g('ssrc')) / ACC, 'SFC LW CRE': (g('str') - g('strc')) / ACC,
              'SFC SW down': g('ssrd') / ACC, 'tcc (model)': g('tcc'), 'lcc (model)': g('lcc')}
    fields['TOA net CRE'] = fields['TOA SW CRE'] + fields['TOA LW CRE']
    fields['SFC net CRE'] = fields['SFC SW CRE'] + fields['SFC LW CRE']
    rows[f'{run} {y0}-{y1}'] = {b: {k: box_mean(v, lat, lon, bx) for k, v in fields.items()} for b, bx in BOXES.items()}

with xr.open_dataset(CER, decode_times=False) as c:
    clat, clon = c['lat'].values, c['lon'].values % 360
    cm = lambda v: c[v].values.mean(0)
    cf = {'TOA SW CRE': cm('toa_cre_sw_clim'), 'TOA LW CRE': cm('toa_cre_lw_clim'),
          'SFC SW CRE': cm('sfc_cre_net_sw_clim'), 'SFC LW CRE': cm('sfc_cre_net_lw_clim'),
          'SFC SW down': cm('sfc_sw_down_all_clim'), 'tcc (model)': cm('cldarea_total_daynight_clim') / 100}
    cf['TOA net CRE'] = cf['TOA SW CRE'] + cf['TOA LW CRE']; cf['SFC net CRE'] = cf['SFC SW CRE'] + cf['SFC LW CRE']
    rows['CERES EBAF (2000s)'] = {b: {k: box_mean(v, clat, clon, bx) for k, v in cf.items()} for b, bx in BOXES.items()}

keys = ['TOA SW CRE', 'TOA LW CRE', 'TOA net CRE', 'SFC SW CRE', 'SFC LW CRE', 'SFC net CRE', 'SFC SW down', 'tcc (model)', 'lcc (model)']
for b in BOXES:
    print(f'\n=== {b}  [W/m2; cloud cover fraction]   (CERES "tcc" row = MODIS cldarea mask, not comparable)')
    print(f'{"":22s}' + ''.join(f'{k:>13s}' for k in keys))
    for name, d in rows.items():
        print(f'{name:22s}' + ''.join(f'{d[b].get(k, np.nan):13.2f}' for k in keys))
