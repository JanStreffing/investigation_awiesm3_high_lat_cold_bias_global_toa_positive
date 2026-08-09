"""Where does the global tropospheric cold bias live? Diagnosis, not a parameter sweep.

WHY.  After 50 runs this is the dominant remaining term for land 2 m temperature: the
boreal-specific excess lives in the bottom three levels, but above 700 hPa Siberia
merely shares a GLOBAL cold bias of 0.7-2.2 K that no lever in the campaign can reach.
It is also the residual in the Siberian budget (the ~0.7 K "universal" component
present over prescribed-SST ocean and under all 20 surface types), the whole of the
2-3.5 degC April-September RIHMI soil deficit, and the prime suspect for the residual
radiative imbalance.  Nobody has asked WHERE it is.

THE QUESTION THAT CHOOSES THE LINE OF ATTACK.  A -2.2 K bias at 300 hPa has very
different implications depending on latitude:

  * TROPICS-LED.  The tropical upper troposphere is tied to the surface by the moist
    adiabat, and SST is PRESCRIBED and correct here.  A cold tropical upper
    troposphere with a correct surface means the convective heating profile is wrong
    -- too much entrainment, too little detrainment aloft, or a convective closure
    problem.  It then propagates globally, because the tropical upper troposphere sets
    the temperature of the whole tropical-to-midlatitude free atmosphere via wave
    dynamics.  That would make this a CONVECTION problem, and it would explain why no
    surface lever has ever touched it.

  * EXTRATROPICS-LED.  Points instead at radiation, the stratospheric boundary
    condition, or resolved dynamics -- a different investigation entirely.

  * UNIFORM.  Points at a radiative or thermodynamic bias acting everywhere: ozone,
    water vapour, or the LW scheme.

THE SECOND QUESTION, from the OLR arithmetic.  toa_decomp_periodclean.py finds the
model 1.6 K cold through the emitting layer but only 2.5 W/m2 short in OLR, where a
Planck scaling predicts about 6.  Roughly 3.6 W/m2 is being offset by something.  A dry
bias would do it -- less water vapour means a lower, warmer effective emission level,
raising OLR.  So q is measured here alongside t, and if the atmosphere is BOTH cold AND
dry, that is a compensating pair: two real errors whose sum is small, which is exactly
what makes a globally tuned budget fragile.

DATA.  amip_presentday 1990-2014 against ERA5 1990-2014, both on the model grid, all 19
model levels present exactly in ERA5's 37 so no vertical interpolation
(vertical_profiles_prep.sh).  Period-clean by construction.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

W = '/tmp/vprof'
LEVELS = [1000, 925, 850, 700, 500, 300, 200, 100]
BANDS = [('90S-60S', -90, -60), ('60S-30S', -60, -30), ('tropics 30S-30N', -30, 30),
         ('30N-60N', 30, 60), ('60N-90N', 60, 90)]

print(__doc__)
print('=' * 100)


def get(src, tag, var, season):
    f = f'{W}/{src}_{tag}_{var}_{season}.nc'
    if not os.path.exists(f):
        return None, None, None
    with xr.open_dataset(f, decode_times=False) as d:
        name = [v for v in d.data_vars if v.lower() in (var, {'t': 'ta', 'q': 'hus',
                                                              'r': 'hur'}.get(var, var))]
        name = name[0] if name else [v for v in d.data_vars if d[v].ndim >= 3][0]
        a = d[name].values
        lat = d['lat'].values
        # the model writes `pressure_levels`, ERA5-via-cdo writes `plev`, and the two
        # are NOT in the same order -- so levels are matched by VALUE below, never by index
        pl = [c for c in ('plev', 'pressure_levels', 'pressure', 'lev', 'level')
              if c in d.coords or c in d.variables]
        lev = np.asarray(d[pl[0]].values, dtype=float) if pl else None
    a = np.squeeze(a)
    if lev is not None and np.nanmax(lev) > 2000:      # Pa -> hPa
        lev = lev / 100.0
    return a, lat, lev


def band(a2d, lat, lo, hi):
    sel = (lat >= lo) & (lat < hi)
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(a2d[sel, :].mean(axis=1), weights=w))


for season in ('ANN', 'JJA', 'DJF'):
    mt, lat, lev = get('model', 'pd', 't', season)
    et, _, elev = get('era5', 'pd', 't', season)
    if mt is None or et is None:
        print(f'{season}: missing input'); continue
    mi = {int(round(p)): i for i, p in enumerate(lev)}
    ei = {int(round(p)): i for i, p in enumerate(elev)}

    print(f'\n{"=" * 100}\n{season}   temperature bias, model - ERA5 [K]\n{"=" * 100}')
    print(f'  {"hPa":>6s} ' + ' '.join(f'{n:>16s}' for n, _, _ in BANDS) + f' {"GLOBAL":>9s}')
    prof = {}
    for p in LEVELS:
        if p not in mi or p not in ei:
            continue
        d = mt[mi[p]] - et[ei[p]]
        vals = [band(d, lat, a, b) for _, a, b in BANDS]
        g = band(d, lat, -90, 90)
        prof[p] = (vals, g)
        print(f'  {p:6d} ' + ' '.join(f'{v:16.2f}' for v in vals) + f' {g:9.2f}')

    if season == 'ANN' and prof:
        print('\n  WHICH BAND CARRIES THE UPPER-TROPOSPHERIC COLD?')
        for p in (300, 200):
            if p not in prof:
                continue
            vals, g = prof[p]
            trop = vals[2]
            others = [vals[0], vals[1], vals[3], vals[4]]
            print(f'    {p} hPa: tropics {trop:+.2f} vs extratropical mean '
                  f'{np.mean(others):+.2f}, global {g:+.2f}')
            spread = max(vals) - min(vals)
            if trop <= min(others):
                v = 'TROPICS-LED -> convective heating profile; SST is prescribed and correct'
            elif trop >= max(others):
                v = 'EXTRATROPICS-LED -> radiation, dynamics or the stratospheric boundary'
            elif spread < 0.6:
                v = 'UNIFORM -> a radiative/thermodynamic bias acting everywhere'
            else:
                v = 'MIXED -> no single band dominates'
            print(f'             band spread {spread:.2f} K  =>  {v}')

# ---------------------------------------------------------------- the dry-bias check
print(f'\n{"=" * 100}')
print('THE COMPENSATING PAIR: is the atmosphere cold AND dry?')
print('=' * 100)
mq, lat, lev = get('model', 'pd', 'q', 'ANN')
eq, _, _ = get('era5', 'pd', 'q', 'ANN')
mt, _, _ = get('model', 'pd', 't', 'ANN')
et, _, _ = get('era5', 'pd', 't', 'ANN')
if mq is not None and eq is not None:
    mi = {int(round(p)): i for i, p in enumerate(lev)}
    _, _, elev2 = get('era5', 'pd', 'q', 'ANN')
    ei = {int(round(p)): i for i, p in enumerate(elev2)}
    print(f'  {"hPa":>6s} {"dT [K]":>9s} {"q model":>10s} {"q ERA5":>10s} '
          f'{"dq [g/kg]":>11s} {"dq / q [%]":>11s}')
    for p in LEVELS:
        if p not in mi or p not in ei:
            continue
        dt = band(mt[mi[p]] - et[ei[p]], lat, -90, 90)
        qm = band(mq[mi[p]], lat, -90, 90) * 1000
        qe = band(eq[ei[p]], lat, -90, 90) * 1000
        rel = 100 * (qm - qe) / qe if qe > 1e-9 else np.nan
        print(f'  {p:6d} {dt:9.2f} {qm:10.4f} {qe:10.4f} {qm - qe:11.4f} {rel:11.1f}')
    # Column water vapour weights the lower troposphere, so a mid/upper dry bias can be
    # radiatively important while barely showing in the total column.
    kk = [p for p in (700, 500, 300) if p in mi and p in ei]
    dq = np.mean([100 * (band(mq[mi[p]], lat, -90, 90) - band(eq[ei[p]], lat, -90, 90))
                  / max(band(eq[ei[p]], lat, -90, 90), 1e-12) for p in kk])
    dT = np.mean([band(mt[mi[p]] - et[ei[p]], lat, -90, 90) for p in kk])
    print(f'\n  700-300 hPa mean:  dT {dT:+.2f} K,  dq {dq:+.1f} % relative')
    if dT < -0.3 and dq < -2:
        print(f'  *** COLD AND DRY.  Both errors are real and they OPPOSE each other in OLR:')
        print(f'      the cold depresses emission, the dryness lowers the effective emission')
        print(f'      level into warmer air and raises it.  That is why a {abs(dT):.1f} K cold')
        print(f'      emitting layer produces only a small OLR deficit, and it is the')
        print(f'      compensating pair the global budget is currently hiding.')
        print(f'      CONSEQUENCE: fixing EITHER one alone moves net TOA the wrong way by')
        print(f'      several W/m2.  They have to be worked together, and any future')
        print(f'      radiation tuning that treats the global mean as nearly-right is')
        print(f'      tuning on top of a cancellation.')
    elif dT < -0.3 and dq > 2:
        print('  Cold and MOIST -- both push OLR down, so the small measured OLR deficit is')
        print('  not explained by water vapour.  Look at cloud LW trapping and cloud-top height.')
    else:
        print('  No clear dry bias in the emitting layer; the OLR residual is elsewhere.')
