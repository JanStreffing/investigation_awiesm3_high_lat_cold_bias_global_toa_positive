"""The missing calibration: Siberian snow cover at GRID scale, from Rutgers 24 km.

WHY THIS EXISTS.  The mode-3 depletion curve was fitted to RIHMI snow courses, which
are 1-2 km transects run only where snow is present.  That constrained the SHAPE of
SCF(d, rho) and the fact that complete cover is real, but it could NOT constrain the
absolute level at a ~100 km grid box, because a course never samples the bare
fraction of a large region.  SCALE was introduced as the free parameter for exactly
that gap and was never calibrated -- P3 uses SCALE=1 (observations taken literally),
P4 SCALE=3, and P3 won on soil temperature.

Rutgers is the dataset that closes it: binary snow/no-snow on a 24 km polar
stereographic grid, so the area fraction covered within the Siberian box is a direct
observational analogue of the model's box-mean cover -- measured over AREA rather
than along a transect.

THE QUESTION IT ANSWERS.  P3 raised Siberian September cover from 0.174 to 0.326 and
October from 0.702 to 0.848 relative to the as-released ramp, while buying only 0.044
of May depletion.  surfece.F90 records that the AS-RELEASED autumn cover was already
+0.18 to +0.25 too high against this same record.  If that still holds, P3 has made a
known autumn error substantially worse in exchange for very little spring, and the
next tuning step is autumn, not spring.

METHOD.  Area-weighted mean of snow_cover_extent over land cells inside
55-75N, 60-180E, by month.  No regridding: the box mask is applied on the native
grid and cells are weighted by their own area, which is what the file provides.

PERIOD CAVEAT.  Rutgers is 1980-2024, the model is pre-industrial.  Autumn onset and
spring melt-out have both shifted over that interval, so the comparison carries a
warming offset of the same kind as every other obs comparison in this campaign.  It
is stated with the numbers rather than corrected for.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

RUT = '/work/ab0246/a270092/obs/snowcover/rutgers_nh_24km_weekly_sce.nc'
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()

# model box-mean cover, from p_series_eval.py (40 yr, same box, same definition)
MODEL = {
    'N1 as-released': [0.964, 0.964, 0.964, 0.944, 0.773, 0.191, None, None,
                       0.174, 0.702, 0.931, 0.962],
    'N2 tanh':        [0.962, 0.963, 0.961, 0.913, 0.615, 0.080, None, None,
                       0.167, 0.680, 0.922, 0.957],
    'P3 fitted':      [0.965, 0.965, 0.964, 0.937, 0.729, 0.167, None, None,
                       0.326, 0.848, 0.957, 0.965],
    'P4 fitted x3':   [0.965, 0.965, 0.963, 0.927, 0.690, 0.149, None, None,
                       0.316, 0.842, 0.958, 0.965],
}

print(__doc__.split('METHOD.')[0])
print('=' * 92)

ds = xr.open_dataset(RUT, decode_times=True)
lat = ds['latitude'].values
lon = ds['longitude'].values
land = ds['land'].values
area = ds['area'].values

# fill values mark cells outside the projection disc
bad = (lat > 90) | (lat < -90) | ~np.isfinite(lat)
box = (~bad) & (lat >= 55) & (lat <= 75) & (((lon + 360) % 360) >= 60) & \
      (((lon + 360) % 360) <= 180) & (land > 0)
print(f'Rutgers cells in the Siberian box over land: {box.sum()}')
print(f'total box area: {area[box].sum():.3e} (file units)\n')

t = ds['time'].values
mo = np.array([int(str(x)[5:7]) for x in t])
w = area[box].astype('float64')

obs = np.full(12, np.nan)
nweeks = np.zeros(12, dtype=int)
for m in range(1, 13):
    k = np.flatnonzero(mo == m)
    if not k.size:
        continue
    tot = 0.0
    for i in k:
        sce = ds['snow_cover_extent'].isel(time=i).values[box].astype('float64')
        sce = np.where(np.isfinite(sce) & (sce <= 1), sce, np.nan)
        tot += np.nansum(sce * w) / np.nansum(np.where(np.isfinite(sce), w, 0.0))
    obs[m - 1] = tot / k.size
    nweeks[m - 1] = k.size
ds.close()

print('OBSERVED Siberian box snow-covered area fraction, Rutgers 24 km 1980-2024\n')
print(f'  {"":16s}' + ''.join(f'{MON[m]:>8s}' for m in range(12)))
print(f'  {"Rutgers":16s}' + ''.join(f'{obs[m]:8.3f}' for m in range(12)))
print(f'  {"(weeks used)":16s}' + ''.join(f'{nweeks[m]:8d}' for m in range(12)))

print('\n\nMODEL - RUTGERS by month  [+ = model holds too much snow]\n')
print(f'  {"run":16s}' + ''.join(f'{MON[m]:>8s}' for m in range(12)))
for lab, v in MODEL.items():
    row = ''.join('     ---' if v[m] is None or not np.isfinite(obs[m])
                  else f'{v[m]-obs[m]:+8.3f}' for m in range(12))
    print(f'  {lab:16s}' + row)

print("""

  READING IT.  The accumulation season (Sep-Oct) and the melt season (May-Jun) are
  separate problems and should be judged separately.

  * If P3's Sep/Oct excess is LARGER than N1's, the fitted curve has worsened the
    autumn error that surfece.F90 already flagged, and the fix is to make d_c larger
    for FRESH low-density snow -- which is precisely what ECE_SNOW_SCF_SWEMIN does.
    SWEMIN was chosen as a numerical stability floor (3.0 kg/m2, the smallest value
    that beats the as-released safety margin); it is ALSO the autumn cover control,
    and that second role was never calibrated against anything.
  * If P3's May/Jun deficit is still positive, the scheme is under-depleting in
    spring and SCALE > 1 is justified after all -- but P4 (SCALE=3) scored WORSE on
    soil temperature, so the two seasons want opposite things and a single SCALE
    cannot serve both.  A larger SWEMIN with a larger SCALE is the combination that
    would fix autumn and spring together, and neither has been tried.""")
