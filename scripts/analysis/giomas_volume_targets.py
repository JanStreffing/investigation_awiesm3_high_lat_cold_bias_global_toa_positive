"""Observed hemispheric sea-ice VOLUME targets, from GIOMAS.

WHY.  The campaign has been scoring Arctic volume against two remembered PIOMAS numbers
(April 28.7, September 13.0) and scoring the Antarctic against nothing at all, because
PIOMAS is Arctic-only.  GIOMAS is the global configuration of the same assimilation
system and carries both hemispheres on one grid, so it gives a consistent target for
each.

WHAT IT IS, AND IS NOT.  GIOMAS is a model with assimilation, not a measurement.  Its
Arctic volume carries roughly 10-20 % uncertainty and its Antarctic volume considerably
more, because Antarctic ice is thin, snow-loaded and poorly constrained by altimetry.
Treat the Antarctic number as an order-of-magnitude anchor, not a target to tune onto.

It is also SATELLITE ERA, so a pre-industrial run is expected to sit ABOVE it, by an
amount this campaign has not pinned down.

heff is effective thickness (volume per unit cell area), so volume = sum(heff * dxt * dyt).

OBS_Y0/OBS_Y1 restrict the window.  The record is 1989-2014 and the Arctic thinned
sharply through it, so the early years are the fairer anchor for a pre-industrial run;
the Antarctic has no significant trend over it and barely moves with the window.

Usage:  python3 scripts/analysis/giomas_volume_targets.py
        OBS_Y0=1989 OBS_Y1=1999 python3 scripts/analysis/giomas_volume_targets.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

F = '/work/ab0246/a270092/obs/GIOMAS/GIOMAS_heff_miss.nc'
with xr.open_dataset(F, decode_times=False) as d:
    heff = np.asarray(d['heff'].values, float)              # (time, j, i) m
    lat  = np.asarray(d['lat_scaler'].values, float)
    dxt  = np.asarray(d['dxt'].values, float)*1e3           # km -> m
    dyt  = np.asarray(d['dyt'].values, float)*1e3
    mon  = np.asarray(d['month'].values, int)
    yr   = np.asarray(d['year'].values, int) if 'year' in d else None
area = dxt*dyt
heff = np.where(np.isfinite(heff), heff, 0.0)
nt = heff.shape[0]
if yr is None:                      # no year coordinate: months run contiguously from 1989
    yr = 1989 + np.arange(nt)//12
Y0 = int(os.environ.get('OBS_Y0', yr.min()))
Y1 = int(os.environ.get('OBS_Y1', yr.max()))
keep = (yr >= Y0) & (yr <= Y1)
print(__doc__.split('\n')[0])
print(f'{int(keep.sum())} months of {nt}, years {Y0}-{Y1}\n')

NH, SH = lat > 0, lat < 0
vol = lambda t, m: float((heff[t]*area*m).sum())/1e12       # m3 -> 10^3 km3

sel = lambda k: [t for t in range(nt) if mon[t] == k and keep[t]]
clim_n = np.array([np.mean([vol(t, NH) for t in sel(k)]) for k in range(1, 13)])
clim_s = np.array([np.mean([vol(t, SH) for t in sel(k)]) for k in range(1, 13)])
sd_n   = np.array([np.std ([vol(t, NH) for t in sel(k)]) for k in range(1, 13)])
sd_s   = np.array([np.std ([vol(t, SH) for t in sel(k)]) for k in range(1, 13)])

names = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
print(f'{"month":<7}{"NH vol":>9}{"+-sd":>7}{"SH vol":>9}{"+-sd":>7}   [10^3 km3]')
for k in range(12):
    print(f'{names[k]:<7}{clim_n[k]:>9.2f}{sd_n[k]:>7.2f}{clim_s[k]:>9.2f}{sd_s[k]:>7.2f}')
print()
print(f'  NH  max {clim_n.max():.2f} ({names[clim_n.argmax()]})   min {clim_n.min():.2f} ({names[clim_n.argmin()]})'
      f'   annual mean {clim_n.mean():.2f}')
print(f'  SH  max {clim_s.max():.2f} ({names[clim_s.argmax()]})   min {clim_s.min():.2f} ({names[clim_s.argmin()]})'
      f'   annual mean {clim_s.mean():.2f}')
