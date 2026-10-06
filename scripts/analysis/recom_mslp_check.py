"""One-month check of the sea-level pressure REcoM uses for gas exchange.

Three runs of January 2140 from the same start:
  PICAL_cc_v350val  tag 2.8.1: REcoM reads press_air, which a coupled FESOM never fills (zero)
  PICAL_cc_slpfix   REcoM uses a constant one atmosphere
  PICAL_cc_mslp     FESOM receives the mean sea-level pressure of OpenIFS into press_air

The pressure REcoM used is not written out, but it can be recovered: dpCO2s is oceanic minus
atmospheric pCO2, and mocsy converts with pCO2atm = xCO2 * (Patm - pH2O).  So
  Patm = (pCO2s - dpCO2s) / xCO2 + pH2O(SST, SSS)
with the Weiss and Price (1980) vapour pressure.  This is compared with the January-mean msl of
OpenIFS (regular grid, nearest neighbour to the FESOM nodes).  Monthly means of products are not
products of monthly means, so agreement is expected to a few tenths of a percent, not exactly.

Usage: recom_mslp_check.py     (reval environment; writes data/clim/recom_mslp_january2140.txt)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import xarray as xr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
RUNS = [('zero', 'PICAL_cc_v350val'), ('constant', 'PICAL_cc_slpfix'), ('coupled', 'PICAL_cc_mslp')]
ATM = 101325.0
BANDS = [('global', -90, 90), ('78-60S', -78, -60), ('60-40S', -60, -40), ('40S-40N', -40, 40),
         ('40-60N', 40, 60), ('60-90N', 60, 90)]

md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lon, lat = md['lon'].values, md['lat'].values
area = md['nod_area'].values[0]


def jan(run, var):
    with xr.open_dataset(f'{P}/runtime/awiesm3-v3.4/{run}/outdata/fesom/{var}.fesom.2140.nc') as ds:
        x = ds[var].sel(time=ds['time'].dt.month == 1).mean('time').values
    return x[:, 0] if x.ndim == 2 and x.shape[1] < x.shape[0] else (x[0] if x.ndim == 2 else x)


def ph2o(t, s):
    tk = t + 273.15
    return np.exp(24.4543 - 67.4509 * (100.0 / tk) - 4.8489 * np.log(tk / 100.0) - 0.000544 * s)


def msl_on_nodes(run):
    with xr.open_dataset(f'{P}/runtime/awiesm3-v3.4/{run}/outdata/oifs/atm_remapped_1m_msl_2140-2140.nc') as ds:
        m = ds['msl'].isel({d: 0 for d in ds['msl'].dims if d not in ('lat', 'lon')})
        return m.sel(lat=xr.DataArray(lat, dims='n'), lon=xr.DataArray(lon % 360, dims='n'), method='nearest').values


def wmean(x, m):
    ok = m & np.isfinite(x)
    return np.sum(x[ok] * area[ok]) / np.sum(area[ok])


D = {}
for tag, run in RUNS:
    d = {v: jan(run, v) for v in ('pCO2s', 'dpCO2s', 'xCO2atm', 'CO2f', 'O2f', 'sst', 'sss')}
    d['pco2atm'] = d['pCO2s'] - d['dpCO2s']
    d['patm'] = d['pco2atm'] / d['xCO2atm'] + ph2o(d['sst'], d['sss'])
    d['msl'] = msl_on_nodes(run) / ATM
    D[tag] = d

wet = np.isfinite(D['coupled']['pCO2s']) & (D['coupled']['xCO2atm'] > 0)
out = ['January 2140, same start. Pressure REcoM used for gas exchange, recovered from pCO2s, dpCO2s, xCO2atm and SST.',
       'zero = PICAL_cc_v350val (tag 2.8.1), constant = PICAL_cc_slpfix (one atmosphere), coupled = PICAL_cc_mslp (OpenIFS msl).',
       '',
       'Recovered pressure [atm], area-weighted, with the OpenIFS January-mean msl of the coupled run:',
       f"{'band':10s} {'zero':>8s} {'constant':>9s} {'coupled':>8s} {'msl':>8s} {'coupled-msl':>12s}"]
for name, a, b in BANDS:
    m = wet & (lat >= a) & (lat < b)
    c = D['coupled']
    out.append(f"{name:10s} {wmean(D['zero']['patm'], m):8.4f} {wmean(D['constant']['patm'], m):9.4f} "
               f"{wmean(c['patm'], m):8.4f} {wmean(c['msl'], m):8.4f} {wmean(c['patm'] - c['msl'], m):+12.4f}")
c = D['coupled']
ok = wet & np.isfinite(c['patm']) & np.isfinite(c['msl'])
r = np.corrcoef(c['patm'][ok], c['msl'][ok])[0, 1]
out += [f"node-by-node, coupled run: correlation of recovered pressure with msl {r:.4f}; "
        f"rms difference {np.sqrt(np.mean((c['patm'][ok] - c['msl'][ok]) ** 2)):.4f} atm; "
        f"range of recovered pressure {np.nanmin(c['patm'][ok]):.4f} to {np.nanmax(c['patm'][ok]):.4f}", '']

out.append('Atmospheric pCO2 seen by the flux [uatm] and what follows from it:')
out.append(f"{'band':10s} {'quantity':28s} {'zero':>9s} {'constant':>9s} {'coupled':>9s} {'cpl-const':>10s}")
for name, a, b in BANDS:
    m = wet & (lat >= a) & (lat < b)
    for key, lab in (('pco2atm', 'atmospheric pCO2 [uatm]'), ('pCO2s', 'surface ocean pCO2 [uatm]'),
                     ('CO2f', 'CO2 flux in [mmolC/m2/d]'), ('O2f', 'O2 flux in [mmolO/m2/d]')):
        v = [wmean(D[t][key], m) for t in ('zero', 'constant', 'coupled')]
        out.append(f"{name:10s} {lab:28s} {v[0]:9.3f} {v[1]:9.3f} {v[2]:9.3f} {v[2] - v[1]:+10.3f}")
out.append('')
for key in ('sst', 'sss'):
    out.append(f"physics, coupled against constant: {key} max |diff| {np.nanmax(np.abs(D['coupled'][key] - D['constant'][key])):.3e}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/recom_mslp_january2140.txt', 'w').write(txt + '\n')
