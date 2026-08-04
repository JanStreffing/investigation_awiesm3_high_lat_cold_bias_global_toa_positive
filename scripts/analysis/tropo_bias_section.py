"""Zonal-mean temperature bias, latitude x pressure: where does the cold troposphere live?

Context. `sub:vprof` established weeks ago that Siberia's free troposphere merely shares a
global tropospheric cold bias, and `vertical_bias_column.py` has now shown that the surface
and soil are the LEAST biased parts of the column while the free troposphere is cold by
~0.9 K and the upper troposphere by ~2.8-3.0 K. With the boreal surface budget nearly spent
(G4 has taken ~0.95 of a 1.3-1.6 K boreal-specific envelope), that tropospheric bias is now
the dominant remaining term for land 2 m temperature -- and it has never been worked on.

This script is the first look. A zonal-mean latitude-pressure section discriminates between
the standard causes far better than any global mean:

  * maximum at the TROPICAL tropopause (~100-200 hPa, 30S-30N)
        -> the classic IFS/EC-Earth cold-point bias: radiation in the UTLS (ozone, water
           vapour), convective detrainment height, or vertical resolution near the cold point.
  * maximum at the POLAR lower stratosphere
        -> polar vortex / gravity-wave drag / ozone.
  * roughly uniform through the troposphere at all latitudes
        -> a bulk radiative or thermodynamic problem, e.g. clear-sky LW or a lapse-rate error.
  * confined to the extratropical lower troposphere
        -> boundary-layer or cloud, i.e. still in reach of the physics this campaign knows.

Land AND ocean together, because the bias exists over prescribed-SST ocean too (-0.72 K at
2 m) and masking to land would hide half the picture.

Period-clean: model `amip_presentday` (1990-2014) vs ERA5 over the same years, ERA5 already
regridded to the model grid by `vertical_bias_column.py` into obs/era5/column/.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

PD = list(range(1990, 2015))
W = '/work/ab0246/a270092/obs/era5/column'
SEAS = {'annual': list(range(12)), 'DJF': [11, 0, 1], 'JJA': [5, 6, 7]}
BANDS = [('90-60S', -90, -60), ('60-30S', -60, -30), ('30S-0', -30, 0),
         ('0-30N', 0, 30), ('30-60N', 30, 60), ('60-90N', 60, 90)]

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def model_pl():
    acc, n, lat, lev = None, 0, None, None
    for y in PD:
        f = f'{RT}/amip_presentday/outdata/oifs/atm_remapped_1m_pl_t_1m_pl_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        d = xr.open_dataset(f)
        a = d['t'].values
        lat = d['t'].lat.values
        lev = d['t'].pressure_levels.values
        d.close()
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, lat, lev) if n else (None, None, None)


def era5_pl():
    p = os.path.join(W, 'e5_t_pl.nc')
    if not os.path.exists(p):
        raise SystemExit('run vertical_bias_column.py first to build the ERA5 column')
    d = xr.open_dataset(p)
    v = [k for k in d.data_vars if d[k].ndim >= 3
         and not np.issubdtype(d[k].dtype, np.datetime64)][0]
    a = d[v].values
    la = d['lat'].values
    lev = d['plev'].values
    d.close()
    if la[0] > la[-1]:
        a = a[:, :, ::-1, :]
    return a, lev


mpl, lat, mlev = model_pl()
if mpl is None:
    raise SystemExit('model pressure-level output missing')
epl, elev = era5_pl()

# match levels by VALUE -- the two files order them differently
pairs = []
for k, p in enumerate(mlev):
    j = int(np.argmin(np.abs(np.asarray(elev) - p)))
    if abs(elev[j] - p) < 1.0:
        pairs.append((p, k, j))
pairs.sort(key=lambda t: -t[0])

w = np.cos(np.deg2rad(lat))

print('Zonal-mean temperature bias, model amip_presentday MINUS ERA5, 1990-2014 [K].')
print('Land AND ocean. Negative = model too cold.\n')
for sname, sidx in SEAS.items():
    print(f'  === {sname} ===')
    print(f'    {"hPa":>7s}' + ''.join(f'{b:>10s}' for b, _, _ in BANDS) + f'{"GLOBAL":>10s}')
    for p, k, j in pairs:
        b = (mpl[sidx, k].mean(0) - epl[sidx, j].mean(0)).mean(axis=1)   # zonal mean
        row = ''
        for _, lo, hi in BANDS:
            m = (lat >= lo) & (lat < hi)
            row += f'{np.average(b[m], weights=w[m]):10.2f}'
        row += f'{np.average(b, weights=w):10.2f}'
        print(f'    {int(round(p/100)):7d}' + row)
    print()

print('  Reading it: a maximum at the tropical tropopause is the classic IFS cold-point')
print('  bias (UTLS radiation / convective detrainment / vertical resolution). A uniform')
print('  tropospheric offset is a bulk radiative or lapse-rate problem. Only a bias')
print('  confined to the extratropical LOWER troposphere is in reach of the surface and')
print('  boundary-layer physics this campaign has been working on.')

# ---------------------------------------------------------------------------
# RESULT (2026-08-04) -- a NEW direction, untouched by all 38 runs.
#
# SHAPE OF THE BIAS. The cold maximises at the TROPOPAUSE at every latitude --
# 200-300 hPa in the extratropics, ~100 hPa in the tropics -- reaching -4.2 K
# (90-60S) and -4.5 K (60-90N) annually, -5.6 K over the Arctic in JJA. The lower
# troposphere is a fairly uniform -0.6 to -1.2 K at all latitudes except 90-60S,
# which is near zero. Above the tropopause the polar sign flips positive.
#
# ENERGETICALLY CONSISTENT, and the cause is NOT cloud:
#     absorbed SW      model 239.98  CERES 241.36   -1.38
#     planetary albedo model 0.2951  CERES 0.2900   +0.0051
#     SW CRE           model -44.42  CERES -45.32   +0.90  <- clouds reflect LESS
#     CLEAR-SKY absorbed SW      284.40 vs 286.68   -2.28  <- the whole shortfall
# Cloud is the wrong sign to explain the deficit. Every lever this campaign built
# was a cloud or surface-vegetation lever; none touched clear-sky shortwave.
#
# AND IT IS AT THE SURFACE, NOT IN THE ATMOSPHERE:
#     TOA clear-sky absorbed        -2.68
#     SURFACE clear-sky absorbed    -3.46
#     atmospheric absorption        +0.78   (wrong sign -- so not aerosol/vapour/ozone)
# Split by surface, per unit area and as a global contribution:
#     land            -6.85 W/m2   (-2.00 global)   albedo +0.030
#     ocean ice-free  -1.38        (-0.89)          albedo +0.005
#     sea ice         -8.59        (-0.57)          albedo +0.056
#
# So the model reflects too much shortwave AT THE SURFACE in clear sky, worst over
# sea ice and land. That is quantitatively consistent with a ~1 K tropospheric cold
# bias and it is a different problem from anything tuned so far.
#
# CAVEATS, and they matter before anyone spends a run:
#  * CERES SURFACE fluxes are a derived product (Kato et al. 2018), constrained by
#    TOA observations plus MODIS surface/cloud properties through a radiative
#    transfer model -- not direct observation. Regional uncertainty is several
#    W/m2, so -6.85 over land is suggestive, not settled.
#  * The albedo columns use the CERES clear-sky DOWNWARD flux as a common
#    denominator for both datasets, because the model does not output clear-sky
#    downward SW at the surface. The measured atmospheric-absorption difference
#    (+0.78 W/m2) makes the resulting error small, but the NET flux differences are
#    the robust numbers and the albedos are indicative.
#  * report sub:albreg found ALL-SKY land albedo essentially perfect over temperate
#    and tropical land (within +-0.004). Reconciling that with +0.030 clear-sky here
#    is the first thing to do -- the two are different quantities and one of them
#    may be wrong.
# ---------------------------------------------------------------------------
