"""CERES or ERA5 -- which is right about Southern Ocean cloud area?

WHY THIS MATTERS ENOUGH TO STOP FOR.  The round-22 headline was that the SO has a
6.43 pp cloud-AMOUNT deficit, and a whole round of runs was designed around closing it.
That number rests on ONE source: CERES cldarea_total_daynight.  Bringing ERA5 in gives
a flatly contradictory answer -- ERA5 SO total cloud is 81.6 % against the model's 83.1,
so against ERA5 the model has slightly TOO MUCH cloud, not too little.  CERES says 89.5.
The two references disagree by 8 pp on the quantity the round was built to fix, and a
single-source inference has already burned this campaign twice (CRUNCEP's 33 W/m2 SW
deficit, the AOD-inferred aerosol term).

WHY A DIRECT COMPARISON OF CLOUD AREAS CANNOT SETTLE IT.  They are not the same
quantity and no amount of care makes them one:

  * CERES cldarea_total_daynight is a MODIS cloud MASK fraction -- the retrieval flags a
    pixel cloudy down to very small optical depth, so optically thin cirrus counts fully.
  * ERA5 tcc and the model's tcc are RADIATIVE total cloud cover, built by overlapping
    layer cloud fractions, which weights by how much cloud is actually there.

A model can therefore be "missing" 6 pp of mask-detectable thin cloud while having
exactly the right radiatively active cover.  Comparing the two numbers is a category
error, and that is what round 22 did.

THE TEST THAT DOES SETTLE IT.  Ask what the cloud DOES rather than how much of it there
is.  Shortwave cloud radiative effect is a difference of two measured broadband fluxes,
identical in definition across model, reanalysis and satellite, and it is the quantity
the tuning actually cares about.  So:

    if ERA5 has less cloud than CERES but reflects the SAME as CERES
        -> the areas differ by DEFINITION, CERES's area is not comparable to a model's,
           and the model's shortfall is in cloud OPACITY, not amount.

    if ERA5 has less cloud AND reflects as little as the model
        -> ERA5 shares the model's bias (both are IFS), CERES is the independent
           arbiter, and the amount deficit is real.

MODIS clt is included as a fourth area estimate.  It is the same instrument family as
the CERES cloud mask, so it is a consistency check on the CERES retrieval rather than an
independent arbiter -- stated so it is not mistaken for one.

UNITS.  ERA5 forecast-stream TOA fields are accumulations.  The divisor is not assumed:
it is derived by requiring ERA5 global absorbed SW to land near 240 W/m2, and the factor
actually used is printed.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

BASE = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive'
W = f'{BASE}/data/vprof'
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
MODISF = '/work/ab0246/a270092/obs/MODIS/clt_MODIS_yearmean.nc'
ACC = 3600.0
SO = (-65.0, -45.0)
BANDS = [('SO 65S-45S', -65, -45), ('90S-65S', -90, -65), ('tropics', -30, 30),
         ('GLOBAL', -90, 90)]

print(__doc__)
print('=' * 96)


def zone(a, lat, lo, hi):
    a = np.squeeze(a)
    if a.ndim == 1:
        a = a[:, None]
    s = (lat >= lo) & (lat < hi)
    return float(np.average(a[s, :].mean(axis=1), weights=np.cos(np.deg2rad(lat[s]))))


def model(var):
    acc, lat = [], None
    for y in range(1990, 2015):
        f = f'{RT}/amip_presentday/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            acc.append(d[var].values.mean(axis=0))
            if lat is None:
                lat = d['lat'].values
    div = 1.0 if var in ('tcc', 'lcc', 'mcc', 'hcc') else ACC
    return np.mean(acc, axis=0) / div, lat


def era(tag, p):
    f = f'{W}/era5_{tag}_{p}.nc'
    if not os.path.exists(f):
        return None, None
    with xr.open_dataset(f, decode_times=False) as d:
        v = [x for x in d.data_vars if d[x].size > 1000][0]
        return np.squeeze(d[v].values), d['lat'].values


cds = xr.open_dataset(CERESF)
clat = cds['lat'].values


def cer(v):
    return cds[v].values.mean(axis=0), clat


# ---- ERA5 accumulation divisor, derived not assumed -----------------------------
e_tsr, elat = era('toa', 178)
if e_tsr is None:
    print('ERA5 TOA files not built yet -- run the cdo prep first.')
    raise SystemExit
raw = zone(e_tsr, elat, -90, 90)
cand = {'already W/m2': 1.0, 'per hour (3600 s)': 3600.0, 'per day (86400 s)': 86400.0}
div = min(cand.items(), key=lambda kv: abs(raw / kv[1] - 240.0))
print(f'ERA5 accumulation: raw global TSR {raw:.4g} -> dividing by {div[0]} gives '
      f'{raw/div[1]:.2f} W/m2 (target ~240)')
D = div[1]

e_ttr, _ = era('toa', 179)
e_tsrc, _ = era('toa', 208)
e_ttrc, _ = era('toa', 209)
e_swcre = (e_tsr - e_tsrc) / D
e_lwcre = (e_ttr - e_ttrc) / D
e_tcc, ectl = era('cloud', 164)

m_tsr, mlat = model('tsr')
m_tsrc, _ = model('tsrc')
m_ttr, _ = model('ttr')
m_ttrc, _ = model('ttrc')
m_tcc, _ = model('tcc')
m_swcre, m_lwcre = m_tsr - m_tsrc, m_ttr - m_ttrc

c_swcre, _ = cer('toa_cre_sw_clim')
c_lwcre, _ = cer('toa_cre_lw_clim')
c_area, _ = cer('cldarea_total_daynight_clim')

with xr.open_dataset(MODISF, decode_times=False) as d:
    mo = np.squeeze(d['clt'].values)
    if mo.ndim == 3:
        mo = mo.mean(axis=0)
    molat = d['lat'].values
if np.nanmax(mo) <= 1.5:
    mo = mo * 100

print(f'\nGlobal sanity: ERA5 SW CRE {zone(e_swcre, elat, -90, 90):.2f} '
      f'(CERES {zone(c_swcre, clat, -90, 90):.2f}), '
      f'ERA5 LW CRE {zone(e_lwcre, elat, -90, 90):.2f} '
      f'(CERES {zone(c_lwcre, clat, -90, 90):.2f})')

print('\n1. CLOUD AREA -- four sources, and they do not measure the same thing')
print('-' * 96)
print(f'  {"band":12s} {"model tcc":>10s} {"ERA5 tcc":>10s} {"CERES mask":>11s} {"MODIS clt":>10s}')
for nm, a, b in BANDS:
    print(f'  {nm:12s} {zone(m_tcc, mlat, a, b)*100:10.1f} {zone(e_tcc, ectl, a, b)*100:10.1f} '
          f'{zone(c_area, clat, a, b):11.1f} {zone(mo, molat, a, b):10.1f}')

print('\n2. WHAT THE CLOUD DOES -- SW cloud radiative effect, identical definition')
print('-' * 96)
print(f'  {"band":12s} {"model":>10s} {"ERA5":>10s} {"CERES":>10s} {"mod-CERES":>11s} '
      f'{"ERA5-CERES":>11s}')
for nm, a, b in BANDS:
    m_, e_, c_ = zone(m_swcre, mlat, a, b), zone(e_swcre, elat, a, b), zone(c_swcre, clat, a, b)
    print(f'  {nm:12s} {m_:10.2f} {e_:10.2f} {c_:10.2f} {m_-c_:+11.2f} {e_-c_:+11.2f}')

print('\n3. THE VERDICT')
print('-' * 96)
so_m = zone(m_swcre, mlat, *SO)
so_e = zone(e_swcre, elat, *SO)
so_c = zone(c_swcre, clat, *SO)
am_m = zone(m_tcc, mlat, *SO) * 100
am_e = zone(e_tcc, ectl, *SO) * 100
am_c = zone(c_area, clat, *SO)
print(f'  Southern Ocean:  area  model {am_m:.1f}  ERA5 {am_e:.1f}  CERES {am_c:.1f} %')
print(f'                   SWCRE model {so_m:.2f}  ERA5 {so_e:.2f}  CERES {so_c:.2f} W/m2')
d_e, d_m = so_e - so_c, so_m - so_c
print(f'  ERA5 reflects {abs(d_e):.2f} W/m2 {"less" if d_e > 0 else "more"} than CERES; '
      f'the model {abs(d_m):.2f} {"less" if d_m > 0 else "more"}.')
print()
if abs(d_e) < 0.4 * abs(d_m):
    print('  *** DEFINITION MISMATCH.  ERA5 has far less cloud AREA than CERES yet reflects')
    print('      almost as much, so CERES cldarea counts optically thin cloud that carries')
    print('      little shortwave.  Its area is NOT comparable to a model tcc, and the')
    print('      "6.43 pp amount deficit" is largely an artefact of that comparison.')
    print('      The model shortfall is real but it is OPACITY: the cloud it has is not')
    print('      bright enough.  That REVERSES the round-23 design, which added cloud.')
elif abs(d_e - d_m) < 0.4 * abs(d_m):
    print('  *** ERA5 SHARES THE MODEL BIAS.  Both are IFS, so this is not independent')
    print('      confirmation of either -- it is the same scheme reproducing itself.')
    print('      CERES is then the only independent arbiter and the amount deficit stands.')
else:
    print('  Mixed: ERA5 sits between CERES and the model, so both an area and an opacity')
    print('  component are present and neither reference can be discarded.')
cds.close()
