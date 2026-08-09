"""Is the extratropical cold-and-dry bias a displaced tropopause?

THE SHAPE IS THE CLUE.  Round 22 found the global cold bias is extratropical, peaks at
200 hPa over both poles (-4.20 Antarctic, -4.45 Arctic) and REVERSES SIGN at 100 hPa
(+0.25, +0.26).  A bias that changes sign across the tropopause is not a bulk cooling --
it is what a DISPLACED tropopause looks like.

  If the model tropopause is too HIGH, then at 200 hPa the model is still on a
  tropospheric lapse rate, cooling with height, while the real atmosphere has already
  turned isothermal.  Model too cold.  Higher up the model finally transitions and the
  difference closes or reverses.  That is the observed pattern.

  If it were too LOW the sign would be the other way round: the model would have warm
  isothermal stratosphere at 200 hPa where reality still has cooling troposphere.

So the hypothesis is testable from the temperature profile alone, without any new run.
This script computes the WMO lapse-rate tropopause for model and ERA5 on identical
levels and compares them per latitude band.

THE SECOND HALF: WHERE THE DRYNESS COMES FROM.  The model is -14.2 % dry at 300 hPa,
-20.2 % at 200 and -46.9 % at 100.  Stratospheric water vapour is set almost entirely by
the temperature of the TROPICAL cold point, which freeze-dries air entering the
stratosphere -- so a cold-biased tropical cold point produces a dry-biased stratosphere
everywhere, including over the poles.  The model's tropical 100 hPa is -2.07 K.  This
script checks whether that is quantitatively enough, via Clausius-Clapeyron on the
saturation mixing ratio at the cold point: if the predicted drying matches the measured
drying, the two halves of the round-22 finding have ONE cause and must be fixed together.

METHOD NOTES.
  * Geometric height increments come from the hypsometric relation, dz = (R Tbar/g)
    ln(p1/p2), computed identically for model and ERA5 so any error is common-mode.
  * WMO definition: lowest level above 500 hPa where the lapse rate falls below 2 K/km.
    The 2 km persistence test cannot be applied faithfully on 19 pressure levels, so it
    is approximated by requiring the NEXT layer to also be below 2 K/km.  That is stated
    rather than hidden -- it makes the absolute tropopause height approximate, but the
    model-minus-ERA5 DIFFERENCE is what the argument uses and both sides get identical
    treatment.
  * Everything is on the model grid, annual mean, amip_presentday 1990-2014 vs ERA5
    1990-2014, so it is period-clean.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

W = os.environ.get('VPROF_DIR', '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/data/vprof')
RD, G = 287.05, 9.80665
BANDS = [('90S-60S', -90, -60), ('60S-30S', -60, -30), ('tropics 30S-30N', -30, 30),
         ('30N-60N', 30, 60), ('60N-90N', 60, 90)]

print(__doc__)
print('=' * 96)


def get(src, var, season='ANN'):
    f = f'{W}/{src}_pd_{var}_{season}.nc'
    with xr.open_dataset(f, decode_times=False) as d:
        name = [v for v in d.data_vars if v.lower() in (var, {'t': 'ta', 'q': 'hus'}.get(var, var))]
        name = name[0] if name else [v for v in d.data_vars if d[v].ndim >= 3][0]
        a = np.squeeze(d[name].values)
        lat = d['lat'].values
        pl = [c for c in ('plev', 'pressure_levels', 'pressure', 'lev', 'level')
              if c in d.coords or c in d.variables][0]
        lev = np.asarray(d[pl].values, dtype=float)
    if np.nanmax(lev) > 2000:
        lev = lev / 100.0
    o = np.argsort(-lev)                      # surface -> top
    return a[o], lat, lev[o]


def bandmean(a2d, lat, lo, hi):
    sel = (lat >= lo) & (lat < hi)
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(a2d[sel, :].mean(axis=1), weights=w))


def profile(a3d, lat, lo, hi):
    return np.array([bandmean(a3d[k], lat, lo, hi) for k in range(a3d.shape[0])])


def tropopause(T, p):
    """WMO lapse-rate tropopause pressure [hPa] from a band-mean profile."""
    lapse = np.full(len(p) - 1, np.nan)
    for k in range(len(p) - 1):
        Tbar = 0.5 * (T[k] + T[k + 1])
        dz = (RD * Tbar / G) * np.log(p[k] / p[k + 1])          # metres, positive up
        lapse[k] = -(T[k + 1] - T[k]) / dz * 1000.0             # K/km, positive = cooling
    for k in range(len(p) - 1):
        if p[k] > 500:                       # WMO: search above 500 hPa only
            continue
        if lapse[k] < 2.0 and (k + 1 >= len(lapse) or lapse[k + 1] < 2.0):
            return p[k], lapse
    return np.nan, lapse


mt, lat, lev = get('model', 't')
et, _, elev = get('era5', 't')
assert np.allclose(lev, elev), 'level mismatch'
mq, _, _ = get('model', 'q')
eq, _, _ = get('era5', 'q')

print(f'levels [hPa]: {", ".join(str(int(round(x))) for x in lev)}\n')
print('1. WMO LAPSE-RATE TROPOPAUSE, model vs ERA5')
print('-' * 96)
print(f'  {"band":18s} {"model":>9s} {"ERA5":>9s} {"diff":>9s}   interpretation')
rows = []
for name, a, b in BANDS:
    Tm, Te = profile(mt, lat, a, b), profile(et, lat, a, b)
    pm, _ = tropopause(Tm, lev)
    pe, _ = tropopause(Te, lev)
    d = pm - pe
    v = ('model tropopause HIGHER (lower p)' if d < -5 else
         'model tropopause LOWER' if d > 5 else 'same to within one level')
    rows.append((name, pm, pe, d))
    print(f'  {name:18s} {pm:9.0f} {pe:9.0f} {d:+9.0f}   {v}')

print('\n2. COLD POINT: temperature and level of the profile minimum')
print('-' * 96)
print(f'  {"band":18s} {"T model":>9s} {"T ERA5":>9s} {"dT":>7s} {"p model":>9s} {"p ERA5":>8s}')
trop_dT = None
for name, a, b in BANDS:
    Tm, Te = profile(mt, lat, a, b), profile(et, lat, a, b)
    im, ie = int(np.argmin(Tm)), int(np.argmin(Te))
    if name.startswith('tropics'):
        trop_dT = Tm[im] - Te[ie]
    print(f'  {name:18s} {Tm[im]:9.2f} {Te[ie]:9.2f} {Tm[im] - Te[ie]:+7.2f} '
          f'{lev[im]:9.0f} {lev[ie]:8.0f}')

# ---------------------------------------------------------------- the drying test
print('\n3. DOES THE TROPICAL COLD POINT EXPLAIN THE STRATOSPHERIC DRYNESS?')
print('-' * 96)


def esat_ice(T):
    """Saturation vapour pressure over ice [Pa], Murphy & Koop (2005)."""
    return np.exp(9.550426 - 5723.265 / T + 3.53068 * np.log(T) - 0.00728332 * T)


Tm_t, Te_t = profile(mt, lat, -30, 30), profile(et, lat, -30, 30)
im, ie = int(np.argmin(Tm_t)), int(np.argmin(Te_t))
pm_cp, pe_cp = lev[im], lev[ie]
qm_pred = 0.622 * esat_ice(Tm_t[im]) / (pm_cp * 100.0)
qe_pred = 0.622 * esat_ice(Te_t[ie]) / (pe_cp * 100.0)
pred_rel = 100 * (qm_pred - qe_pred) / qe_pred
print(f'  tropical cold point: model {Tm_t[im]:.2f} K at {pm_cp:.0f} hPa, '
      f'ERA5 {Te_t[ie]:.2f} K at {pe_cp:.0f} hPa  ({trop_dT:+.2f} K)')
print(f'  saturation mixing ratio there: model {qm_pred*1e6:.3f} vs ERA5 {qe_pred*1e6:.3f} ppmv-mass')
print(f'  => Clausius-Clapeyron predicts entry-air drying of {pred_rel:+.1f} %')

k100 = int(np.argmin(np.abs(lev - 100)))
meas = 100 * (bandmean(mq[k100], lat, -90, 90) - bandmean(eq[k100], lat, -90, 90)) \
       / bandmean(eq[k100], lat, -90, 90)
print(f'  measured global drying at 100 hPa: {meas:+.1f} %')
print()
if pred_rel < 0 and meas < 0:
    frac = pred_rel / meas
    print(f'  The cold point accounts for {100*frac:.0f} % of the measured drying.')
    if 0.5 <= frac <= 1.8:
        print('  *** ONE CAUSE.  A cold-biased tropical cold point freeze-dries the air')
        print('      entering the stratosphere, and that dry air is transported poleward by')
        print('      the Brewer-Dobson circulation -- so the polar UTLS dryness is not a')
        print('      separate error, it is the tropical cold point seen downstream.')
        print('      CONSEQUENCE: cold and dry are ONE problem with ONE entry point, and a')
        print('      lever that warms the tropical cold point should fix both.  That is a')
        print('      much better prospect than two independent errors that cancel.')
    else:
        print('      Only partial: the cold point contributes but does not close it, so')
        print('      there is a second dry source (transport, or the model\'s own UTLS')
        print('      microphysics) still to find.')
else:
    print('  Signs do not line up; the cold point is not the dryness source.')

print('\n4. LAPSE RATE THROUGH THE UTLS, model - ERA5 [K/km]  (positive = model cools faster)')
print('-' * 96)
hdr = [f'{int(round(lev[k]))}-{int(round(lev[k+1]))}' for k in range(len(lev) - 1)
       if 400 >= lev[k + 1] >= 50]
print(f'  {"band":18s} ' + ' '.join(f'{h:>10s}' for h in hdr))
for name, a, b in BANDS:
    Tm, Te = profile(mt, lat, a, b), profile(et, lat, a, b)
    _, lm = tropopause(Tm, lev)
    _, le = tropopause(Te, lev)
    vals = [lm[k] - le[k] for k in range(len(lev) - 1) if 400 >= lev[k + 1] >= 50]
    print(f'  {name:18s} ' + ' '.join(f'{v:10.2f}' for v in vals))
print('\n  A model that cools faster than ERA5 through the layer BELOW its tropopause and')
print('  slower above it is the quantitative statement of "tropopause too high".')
