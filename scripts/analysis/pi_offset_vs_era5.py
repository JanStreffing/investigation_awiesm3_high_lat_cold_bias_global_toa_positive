"""Is the run cold-biased, or is it just pre-industrial?

Every observational yardstick this campaign uses is a present-day one: ERA5 T2m 1990-2014,
CERES 2000-2020, HadISST satellite era.  A PI-forced run is SUPPOSED to sit below all of
them, so a raw "bias" against any of them is a bias plus the industrial-era warming, and
reading it as error alone manufactures a cold bias that is not there.

AR6 puts GSAT 1850-1900 -> 1995-2014 at +0.85 K [0.67-0.98], and the ERA5 reference here is
1990-2014, so the expected global offset for a correct PI run is about -0.85 K.  Arctic
amplification of 2-4x makes the expected 60-90N offset roughly -1.7 to -3.4 K annual, and
larger still in DJF because the Arctic warming is winter-concentrated.

This prints the model minus ERA5 offsets by band and season so they can be read against
those expectations rather than against zero.  It does NOT tell you the model is right --
the expected offsets carry real uncertainty, and the Arctic DJF one especially.  It tells
you how much of an apparent bias is accounted for before any error is claimed.

Usage:  Y0=1950 Y1=1959 ARM=PICAL_momixoff python3 scripts/analysis/pi_offset_vs_era5.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob, warnings
import numpy as np, xarray as xr
warnings.filterwarnings('ignore')

R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
E5T2M = '/work/ab0246/a270092/obs/era5/netcdf/T2M.nc'   # 1990-2014
ARM = os.environ.get('ARM', 'PICAL_momixoff')
Y0, Y1 = int(os.environ.get('Y0', 1950)), int(os.environ.get('Y1', 1959))
SEAS = {'ANN': list(range(12)), 'DJF': [11, 0, 1], 'JJA': [5, 6, 7]}

# AMIP arms write no lsm; fall back to the coupled grid, which is the same TCO95 mesh.
_lsm = (sorted(glob.glob(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_lsm_*.nc'))
        or sorted(glob.glob('/work/bb1469/a270092/runtime/awiesm3-v3.4/PICAL_momixoff/'
                            'outdata/oifs/atm_remapped_1m_lsm_*.nc')))
f0 = _lsm[0]
with xr.open_dataset(f0, decode_times=False) as d:
    msk = np.squeeze(d['lsm'].values); msk = msk[0] if msk.ndim == 3 else msk
    lat = np.squeeze(d['lat'].values); lon = np.squeeze(d['lon'].values)
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], msk.shape)
LAT = np.broadcast_to(lat[:, None], msk.shape)

acc = []
for y in range(Y0, Y1 + 1):
    # coupled arms write atm_remapped_1m_2t_YYYY-YYYY.nc, AMIP arms ..._2t_1m_YYYY-YYYY.nc
    g = (glob.glob(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc')
         or glob.glob(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_2t_1m_{y}-{y}.nc'))
    with xr.open_dataset(g[0], decode_times=False) as d:
        k2 = '2t' if '2t' in d.data_vars else [c for c in d.data_vars if 'bnds' not in c][0]
        acc.append(np.squeeze(d[k2].values))
mod = np.mean(acc, axis=0)                                   # (12, lat, lon)
mod = mod - 273.15 if np.nanmean(mod) > 100 else mod

with xr.open_dataset(E5T2M, decode_times=True) as d:
    clim = d['t2m'].groupby('time.month').mean('time')
    la = [c for c in clim.dims if 'lat' in c][0]; lo = [c for c in clim.dims if 'lon' in c][0]
    obs = np.asarray(clim.interp({la: ('lat', lat), lo: ('lon', lon)}).values, dtype=float)
obs = obs - 273.15 if np.nanmean(obs) > 100 else obs
nh = lat > 30                                                # month-axis assertion
cyc = [float(np.average(obs[i][nh].mean(1), weights=np.cos(np.deg2rad(lat))[nh])) for i in range(12)]
if int(np.argmax(cyc)) not in (5, 6, 7):
    raise SystemExit(f'ERA5 month axis wrong: NH peak at index {np.argmax(cyc)}')

BANDS = [('GLOBAL', -90, 90, None), ('60-90N', 60, 90, None), ('60-90N land', 60, 90, True),
         ('Siberia 55-75N 60-180E', 55, 75, True), ('30-60N', 30, 60, None),
         ('20S-20N', -20, 20, None), ('60-90S', -90, -60, None)]
EXPECT = {'GLOBAL': '-0.85 (AR6 GSAT 1850-1900 -> 1995-2014)',
          '60-90N': '-1.7 to -3.4 annual (2-4x amplification); more in DJF',
          '60-90N land': 'as 60-90N', '30-60N': '-0.6 to -1.2', '20S-20N': '-0.6 to -0.9',
          '60-90S': '-0.5 to -1.0 (weak southern amplification)',
          'Siberia 55-75N 60-180E': '-1.5 to -2.5 annual, -1.0 to -1.5 JJA (soft: my\n                                    construction from regional rates, not a cited value)'}

print(__doc__.split('Usage:')[0])
print(f'{ARM}  {Y0}-{Y1}   model minus ERA5 1990-2014 [K]\n'
      '(for a PRESENT-DAY-forced run the expected offset column is void: read against zero)\n')
print(f'{"band":<24}{"ANN":>8}{"DJF":>8}{"JJA":>8}   expected PI offset')
for name, la_, lb_, land in BANDS:
    k = (LAT >= la_) & (LAT <= lb_)
    if land: k = k & (msk > 0.5)
    if 'Siberia' in name: k = k & (((lon[None, :] - 60) % 360) <= 120)   # canonical campaign box
    row = []
    for s in ('ANN', 'DJF', 'JJA'):
        d = (mod[SEAS[s]].mean(0) - obs[SEAS[s]].mean(0))
        row.append(float(np.average(d[k], weights=W[k])))
    print(f'{name:<24}{row[0]:8.2f}{row[1]:8.2f}{row[2]:8.2f}   {EXPECT[name]}')
print('\nA band sitting AT its expected offset carries no detectable bias; one well below it\n'
      'is genuinely cold.  The expectations are literature ranges, not measurements.')
