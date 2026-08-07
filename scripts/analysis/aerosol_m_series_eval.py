"""M series: how much of the clear-sky reflection excess is aerosol?

WHY THIS EXISTS.  The campaign's energy problem decomposes into TWO reflection
columns, and they must never be collapsed into the all-sky total because in the
Southern Ocean they nearly cancel (+5.10 clear-sky against -7.65 cloud gives only
-2.54 all-sky, so the SO looks nearly right while both halves are badly wrong):

                      clear-sky   cloud
    60-90N              +4.42     +1.36
    30-60N              +4.22     -4.22
    tropics             +1.55     +0.41
    SO 45-65S           +5.10     -7.65
    60-90S              +2.34     -8.19
    GLOBAL              +2.68     -1.17      [W/m2 vs CERES]

Of the SO +5.10, roughly 2.2 was attributed to sea-ice albedo, 0.9 to ocean surface
albedo, ~0.7 to aerosol (from a +0.033 AOD excess vs MISR) and ~1.3 left unexplained.
The aerosol term was an INFERENCE from an AOD comparison, never a model experiment.
M1 and M2 test it directly.

    M1  LMACV2SP = .false.   anthropogenic MACv2-SP plumes removed entirely
    M2  LAER3D   = .true.    CAMS climatology given its 3D vertical distribution
                             instead of the 2D column

BOTH ARE ON THE PRESENT-DAY BASE (1989-2014), not the PI base.  MACv2-SP is
transient, so at 1872-1915 the anthropogenic plumes are already near zero and the
test would show nothing.  Present-day is also the period the +2.68 W/m2 was measured
over and the only one comparable with CERES.  They are in NOT_LEVERS so evaluate.sh
cannot difference them against the PI control over the wrong years; the comparison
here is against amip_presentday over the SAME years.

THE DECOMPOSITION.  With tisr the incident SW, tsr the all-sky net and tsrc the
clear-sky net (all TOA, positive down):
    clear-sky reflection = tisr - tsrc
    cloud reflection     = tsrc - tsr        (= -CRE_SW)
    all-sky reflection   = tisr - tsr        (the sum -- reported but never used alone)

WHAT TO READ
  * M1 - presentday on the CLEAR-SKY column is the anthropogenic aerosol's direct
    contribution, with sign flipped: removing the aerosol should REDUCE clear-sky
    reflection, so a negative delta of magnitude X means aerosol contributes +X.
    Compare that against the ~0.7 W/m2 the AOD inference predicted for the SO.
  * M1 on the CLOUD column is the indirect (aerosol-cloud) effect, which
    LMACV2SP_CCNF routes through the CCN field.  If it is large, the two columns
    are not separable and the attribution table needs rewriting.
  * M2 - presentday isolates vertical redistribution at fixed mass.  The 2D and 3D
    CAMS files carry the same column mass to 0.4 %, so any signal here is purely
    where the aerosol sits, not how much there is.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

from runs import RT, OBS

# IFS TOA fluxes are ACCUMULATED over the output step, in J/m2, not W/m2.  With a
# 3600 s timestep the raw numbers come out ~3600x too large (global clear-sky
# reflection 203270 instead of 56).  Verified against the campaign's own documented
# table: /3600 gives SO clear-sky +4.8 and cloud -7.5 vs the recorded +5.10/-7.65.
ACC = 3600.0

YEARS = range(1990, 2015)          # period-clean overlap with CERES-era forcing
RUNS = [('presentday', 'amip_presentday'),
        ('M1 no-anth', 'amip_M1_noanthaer'),
        ('M2 aer3d',   'amip_M2_aer3d')]
BANDS = [('60-90N', 60, 90), ('30-60N', 30, 60), ('tropics', -30, 30),
         ('SO 45-65S', -65, -45), ('60-90S', -90, -60), ('GLOBAL', -90, 90)]


def monthly(run, var):
    acc, n = None, 0
    for y in YEARS:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values
            la, lo = d['lat'].values, d['lon'].values
        if a.shape[0] != 12:
            continue
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, la, lo, n) if n else (None, None, None, 0)


print(__doc__.split('WHAT TO READ')[0])
print('=' * 100)

# ---- model: build both reflection columns -------------------------------------
M = {}
for lab, run in RUNS:
    tisr, la, lo, n = monthly(run, 'tisr')
    tsr, _, _, _ = monthly(run, 'tsr')
    tsrc, _, _, _ = monthly(run, 'tsrc')
    if tisr is None or tsr is None or tsrc is None:
        print(f'  {lab}: missing TOA output'); continue
    M[lab] = dict(clr=(tisr - tsrc)/ACC, cld=(tsrc - tsr)/ACC,
                  all=(tisr - tsr)/ACC, lat=la, n=n)
    print(f'  {lab:12s} {n} yr')

LAT, LON = M['presentday']['lat'], M['presentday']['lon'] if 'lon' in M['presentday'] else None
w2 = np.cos(np.deg2rad(LAT))


def band(field, lo_, hi_):
    """Annual, area-weighted mean over a latitude band."""
    k = (LAT >= lo_) & (LAT <= hi_)
    a = field.mean(axis=0)                       # annual mean, (nlat, nlon)
    ww = np.broadcast_to(w2[:, None], a.shape)[k]
    return float(np.average(a[k], weights=ww))


# ---- CERES on the same decomposition ------------------------------------------
with xr.open_dataset(OBS) as d:
    csw = d['toa_sw_all_clim'].values            # all-sky reflected SW (up)
    csw_clr = d['toa_sw_clr_t_clim'].values      # clear-sky reflected, total-region
    clat = d['lat'].values
cw = np.cos(np.deg2rad(clat))


def cband(field, lo_, hi_):
    k = (clat >= lo_) & (clat <= hi_)
    a = field.mean(axis=0)
    ww = np.broadcast_to(cw[:, None], a.shape)[k]
    return float(np.average(a[k], weights=ww))


print('\n\n1. ABSOLUTE reflection by band [W/m2]  (clear-sky | cloud), model vs CERES\n')
hdr = f'  {"band":11s}' + ''.join(f'{l:>21s}' for l, _ in [(x, 0) for x in
      ('CERES', 'presentday', 'M1 no-anth', 'M2 aer3d')])
print(f'  {"band":11s}' + ''.join(f'{n:>22s}' for n in ('CERES', 'presentday', 'M1 no-anth', 'M2 aer3d')))
print(f'  {"":11s}' + ''.join(f'{"clr":>11s}{"cld":>11s}' for _ in range(4)))
for name, lo_, hi_ in BANDS:
    row = f'  {name:11s}'
    ccl = cband(csw_clr, lo_, hi_); cal = cband(csw, lo_, hi_)
    row += f'{ccl:11.2f}{cal-ccl:11.2f}'
    for lab in ('presentday', 'M1 no-anth', 'M2 aer3d'):
        if lab not in M:
            row += f'{"":11s}{"":11s}'; continue
        row += f'{band(M[lab]["clr"], lo_, hi_):11.2f}{band(M[lab]["cld"], lo_, hi_):11.2f}'
    print(row)

print('\n\n2. MODEL - CERES bias by band [W/m2]   (this is the +2.68 / -1.17 table)\n')
print(f'  {"band":11s}' + ''.join(f'{n:>24s}' for n in ('presentday', 'M1 no-anth', 'M2 aer3d')))
print(f'  {"":11s}' + ''.join(f'{"clr":>12s}{"cld":>12s}' for _ in range(3)))
for name, lo_, hi_ in BANDS:
    ccl = cband(csw_clr, lo_, hi_); ccd = cband(csw, lo_, hi_) - ccl
    row = f'  {name:11s}'
    for lab in ('presentday', 'M1 no-anth', 'M2 aer3d'):
        if lab not in M:
            row += f'{"":12s}{"":12s}'; continue
        row += f'{band(M[lab]["clr"], lo_, hi_)-ccl:+12.2f}{band(M[lab]["cld"], lo_, hi_)-ccd:+12.2f}'
    print(row)

print('\n\n3. THE EXPERIMENT: M - presentday, both columns [W/m2]\n')
print('   M1 clear-sky delta = MINUS the anthropogenic aerosol direct effect.')
print('   AOD inference predicted aerosol contributes ~0.7 W/m2 of the SO +5.10.\n')
print(f'  {"band":11s}' + ''.join(f'{n:>26s}' for n in ('M1 - presentday', 'M2 - presentday')))
print(f'  {"":11s}' + ''.join(f'{"clr":>13s}{"cld":>13s}' for _ in range(2)))
for name, lo_, hi_ in BANDS:
    row = f'  {name:11s}'
    for lab in ('M1 no-anth', 'M2 aer3d'):
        if lab not in M:
            row += f'{"":13s}{"":13s}'; continue
        row += (f'{band(M[lab]["clr"],lo_,hi_)-band(M["presentday"]["clr"],lo_,hi_):+13.3f}'
                f'{band(M[lab]["cld"],lo_,hi_)-band(M["presentday"]["cld"],lo_,hi_):+13.3f}')
    print(row)

print("""

  READING IT.  In section 3, a NEGATIVE M1 clear-sky delta means removing the
  anthropogenic aerosol removed that much reflection, i.e. the aerosol was
  contributing +|delta| to the clear-sky excess.  If that number is ~0.7 in the SO
  the AOD inference is confirmed; if it is much smaller, the aerosol term in the
  attribution table is wrong and the ~1.3 unexplained residual is really ~2.0.
  A large M1 CLOUD delta means the direct and indirect effects are not separable
  and the two-column attribution needs rebuilding.  M2 moves mass in the vertical
  only, so a large signal there would mean the 2D CAMS assumption -- not the
  aerosol amount -- is a first-order error.""")
