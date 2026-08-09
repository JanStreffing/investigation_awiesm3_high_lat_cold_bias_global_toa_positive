"""Anatomy of the Southern Ocean cloud deficit -- the largest unexploited term left.

WHY.  The SO band (65-45S) absorbs +7.5 W/m2 too much energy period-cleanly
(toa_decomp_periodclean.py), contributing ~+0.7 W/m2 to the global imbalance on its
own, and the cause is a cloud deficit: area 6.4 pp below CERES with SW CRE 7.5 W/m2
too weak.  Aerosol has been measured OUT of it (-0.116 W/m2, six times smaller than
the AOD inference), so the residual is now LARGER than when the campaign started
looking.  Every lever so far has moved it only as a by-product.

Before designing a round, four things have to be known, and none of them can be read
off a CRE number:

  1. AMOUNT OR OPACITY.  CRE = area x (CRE per unit area).  A cloud-amount lever and a
     cloud-optics lever are different physics with different side effects, and the
     campaign has been scoring both on the same aggregate.  Splitting them tells you
     which knob can even reach the error.  The split is done at OBSERVED CRE per unit
     area, so the "amount" term is what the missing cloud would contribute if it
     looked like the cloud that is there.

  2. WHERE.  Zonally uniform, or concentrated in a sector?  A sector-localised deficit
     points at a circulation or sea-ice-edge problem; a uniform one points at cloud
     microphysics.  The D2b result (ice nuclei) predicts uniform.

  3. WHEN.  The supercooled-liquid hypothesis predicts the deficit peaks in austral
     summer, when the SW is there to be reflected and mixed-phase cloud is most
     exposed; a wintertime-dominated deficit would mean something else.

  4. WHAT THE LEVERS ACTUALLY DID TO AREA.  D2a/D2b were adopted on SW RMSE and CRE.
     If they closed CRE by making existing cloud brighter rather than by making more
     of it, then the area deficit is untouched and the mechanism is not what we think.
     This is the falsifiable part, and it is scored here on AREA for the first time.

VERTICAL PLACEMENT is included where the model supports it: lcc/mcc/hcc exist in the
monthly output, and CERES supplies cloud effective pressure, so "is the model's SO
cloud in the right layer" is answerable even though the model writes no cloud fraction
on pressure levels.

PERIOD.  amip_presentday 1990-2014 against CERES 07/2005-06/2015 for the control
anatomy -- period-clean, and the epoch offset on the SO cloud term was measured at
-0.07 W/m2, i.e. negligible.  The LEVER comparison necessarily uses the PI-epoch
arms, because that is where D2a/D2b were run; lever DELTAS are epoch-insensitive since
both sides share the forcing.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

ACC = 3600.0
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
SO = (-65.0, -45.0)
PD, PDY = 'amip_presentday', (1990, 2014)
# lever arms: all PI-epoch, all 44 yr, deltas against the shared control
# EVERY run that set a NAMCLDP parameter, so the "has this been tried" question is
# answered by the table instead of from memory.  Getting that wrong cost a duplicate
# run on 2026-08-09: RCL_OVERLAPLIQICE=0.35 was proposed as new when it is A1b, and
# "ovl" in the A-series labels means the liquid/ice DEPOSITION overlap, not the
# radiative cloud overlap.  The two land/sea levers at the bottom are the control:
# they cannot physically act on Southern Ocean cloud, so whatever they register is
# the noise floor of this diagnostic.
LEVERS = [('control', 'amip_pi_base'),
          ('A1a ovlliqice=0.10', 'amip_A1_overlap01'),
          ('A1b ovlliqice=0.35', 'amip_A1_overlap035'),
          ('A1c depliqdep1500', 'amip_A1c_depliqdepth1500'),
          ('B2 clddiff_convi25', 'amip_B2_clddiffconvi25'),
          ('B3 clddiff 1.5e-5', 'amip_B3_clddiff15e6'),
          ('B6 lcritsnow 1e-5', 'amip_B6_lcritsnow1e5'),
          ('B7 rvice 0.22', 'amip_B7_rvice022'),
          ('D2a inpsea 0.2', 'amip_D2a_inpsea02'),
          ('D2b inp+p700', 'amip_D2b_inpsea02_p700'),
          ('G4 tundra225 [land]', 'amip_G4_tundra'),
          ('K1 landalb [land]', 'amip_K1_landalb')]
LEVY = (1872, 1915)
SEAS = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8], 'SON': [9, 10, 11]}

print(__doc__)
print('=' * 100)


# ------------------------------------------------------------------ model loading
def monthly(run, var, y0, y1):
    """(12, nlat, nlon) monthly climatology, and lat.  Radiative vars are de-accumulated."""
    D = f'{RT}/{run}/outdata/oifs'
    acc, lat = [], None
    div = 1.0 if var in ('tcc', 'lcc', 'mcc', 'hcc') else ACC
    for y in range(y0, y1 + 1):
        f = f'{D}/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values / div
            if lat is None:
                lat = d['lat'].values
        if a.shape[0] == 12:
            acc.append(a)
    if not acc:
        return None, None, 0
    return np.mean(acc, axis=0), lat, len(acc)


def zone(field2d, lat, a=SO[0], b=SO[1]):
    sel = (lat >= a) & (lat < b)
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(field2d[sel, :].mean(axis=1), weights=w))


# ------------------------------------------------------------------ CERES
cds = xr.open_dataset(CERESF)
clat, clon = cds['lat'].values, cds['lon'].values
csel = (clat >= SO[0]) & (clat < SO[1])
cw = np.cos(np.deg2rad(clat[csel]))


def cmon(v):
    return cds[v].values          # (12, lat, lon)


def czone(v, months=None):
    a = cmon(v)
    if months is not None:
        a = a[[m - 1 for m in months]]
    return float(np.average(a.mean(axis=0)[csel, :].mean(axis=1), weights=cw))


def czone_lon(v, months=None):
    """SO-band mean as a function of longitude."""
    a = cmon(v)
    if months is not None:
        a = a[[m - 1 for m in months]]
    return clon, np.average(a.mean(axis=0)[csel, :], axis=0, weights=cw)


# ------------------------------------------------------------------ control anatomy
tcc, lat, ny = monthly(PD, 'tcc', *PDY)
tsr, _, _ = monthly(PD, 'tsr', *PDY)
tsrc, _, _ = monthly(PD, 'tsrc', *PDY)
ttr, _, _ = monthly(PD, 'ttr', *PDY)
ttrc, _, _ = monthly(PD, 'ttrc', *PDY)
lcc, _, _ = monthly(PD, 'lcc', *PDY)
mcc, _, _ = monthly(PD, 'mcc', *PDY)
hcc, _, _ = monthly(PD, 'hcc', *PDY)
swcre = tsr - tsrc
print(f'model {PD} {PDY[0]}-{PDY[1]}: {ny} years; SO band {SO[0]:.0f} to {SO[1]:.0f}\n')

m_area = zone(tcc.mean(axis=0), lat) * 100
o_area = czone('cldarea_total_daynight_clim')
m_cre = zone(swcre.mean(axis=0), lat)
o_cre = czone('toa_cre_sw_clim')

print('1. AMOUNT OR OPACITY -- splitting the SW CRE error')
print('-' * 100)
print(f'  {"":26s} {"model":>9s} {"CERES":>9s} {"diff":>9s}')
print(f'  {"cloud area [%]":26s} {m_area:9.2f} {o_area:9.2f} {m_area - o_area:+9.2f}')
print(f'  {"SW CRE [W/m2]":26s} {m_cre:9.2f} {o_cre:9.2f} {m_cre - o_cre:+9.2f}')
cre_per_area_o = o_cre / (o_area / 100)
cre_per_area_m = m_cre / (m_area / 100)
print(f'  {"CRE per unit area":26s} {cre_per_area_m:9.2f} {cre_per_area_o:9.2f} '
      f'{cre_per_area_m - cre_per_area_o:+9.2f}')
amount_term = (m_area - o_area) / 100 * cre_per_area_o
opacity_term = (m_cre - o_cre) - amount_term
print(f'\n  decomposition of the {m_cre - o_cre:+.2f} W/m2 CRE error:')
print(f'    AMOUNT  (missing area x observed CRE/area) {amount_term:+8.2f} W/m2  '
      f'{100 * amount_term / (m_cre - o_cre):5.1f} %')
print(f'    OPACITY (residual)                         {opacity_term:+8.2f} W/m2  '
      f'{100 * opacity_term / (m_cre - o_cre):5.1f} %')
dom = 'AMOUNT' if abs(amount_term) > abs(opacity_term) else 'OPACITY'
print(f'  -> {dom}-dominated.  ' + ('A cloud-amount lever can reach most of this.'
      if dom == 'AMOUNT' else 'A cloud-amount lever alone CANNOT reach most of this.'))

print('\n2. WHERE -- SO cloud-area error by longitude sector')
print('-' * 100)
_, o_lon = czone_lon('cldarea_total_daynight_clim')
mlon = None
with xr.open_dataset(f'{RT}/{PD}/outdata/oifs/atm_remapped_1m_tcc_1m_2000-2000.nc',
                     decode_times=False) as d:
    mlon = d['lon'].values
msel = (lat >= SO[0]) & (lat < SO[1])
mw = np.cos(np.deg2rad(lat[msel]))
m_lon = np.average(tcc.mean(axis=0)[msel, :], axis=0, weights=mw) * 100
SECT = [('Atlantic  60W-20E', -60, 20), ('Indian    20E-90E', 20, 90),
        ('Australia 90E-180', 90, 180), ('Pacific   180-60W', 180, 300)]
print(f'  {"sector":22s} {"model":>8s} {"CERES":>8s} {"diff":>8s}')
for name, a, b in SECT:
    mm = ((mlon - a) % 360) < ((b - a) % 360 or 360)
    cc = ((clon - a) % 360) < ((b - a) % 360 or 360)
    print(f'  {name:22s} {m_lon[mm].mean():8.2f} {o_lon[cc].mean():8.2f} '
          f'{m_lon[mm].mean() - o_lon[cc].mean():+8.2f}')
spread = np.ptp([m_lon[((mlon - a) % 360) < ((b - a) % 360 or 360)].mean()
                 - o_lon[((clon - a) % 360) < ((b - a) % 360 or 360)].mean()
                 for _, a, b in SECT])
print(f'  sector spread {spread:.2f} pp -> ' +
      ('ZONALLY UNIFORM: consistent with microphysics, not circulation.' if spread < 3
       else 'SECTOR-DEPENDENT: circulation or ice edge is involved, not microphysics alone.'))

print('\n3. WHEN -- by season, and against the shortwave that is available to reflect')
print('-' * 100)
print(f'  {"season":8s} {"m area":>8s} {"o area":>8s} {"d area":>8s} {"m CRE":>8s} '
      f'{"o CRE":>8s} {"d CRE":>8s} {"d net":>8s}')
for s, mo in SEAS.items():
    idx = [m - 1 for m in mo]
    ma = zone(tcc[idx].mean(axis=0), lat) * 100
    oa = czone('cldarea_total_daynight_clim', mo)
    mc = zone(swcre[idx].mean(axis=0), lat)
    oc = czone('toa_cre_sw_clim', mo)
    mn = zone((tsr + ttr)[idx].mean(axis=0), lat)
    on = czone('toa_net_all_clim', mo)
    print(f'  {s:8s} {ma:8.2f} {oa:8.2f} {ma - oa:+8.2f} {mc:8.2f} {oc:8.2f} '
          f'{mc - oc:+8.2f} {mn - on:+8.2f}')

print('\n4. VERTICAL PLACEMENT -- which layer holds the SO cloud')
print('-' * 100)
print(f'  {"model low  (lcc)":26s} {zone(lcc.mean(axis=0), lat) * 100:9.2f} %')
print(f'  {"model mid  (mcc)":26s} {zone(mcc.mean(axis=0), lat) * 100:9.2f} %')
print(f'  {"model high (hcc)":26s} {zone(hcc.mean(axis=0), lat) * 100:9.2f} %')
print(f'  {"model total (tcc)":26s} {m_area:9.2f} %')
mp = czone('cldpress_total_daynight_clim')
print(f'  {"CERES eff. cloud pressure":26s} {mp:9.1f} hPa   '
      f'(low cloud if >680, mid 440-680, high <440)')
print(f'  {"CERES cloud optical depth":26s} {czone("cldtau_total_day_clim"):9.2f}')
print('  The D2b pressure gate acts below 700 hPa, i.e. on the LOW branch; if CERES puts')
print('  the effective cloud there too, the gate is aimed at the right layer.')

# ------------------------------------------------------------------ 5. levers on AREA
print('\n5. WHAT THE LEVERS DID TO CLOUD AREA (not CRE) -- the falsifiable part')
print('-' * 100)
base = None
rows = []
for name, run in LEVERS:
    t, la, n = monthly(run, 'tcc', *LEVY)
    if t is None:
        print(f'  {name:16s} no tcc output'); continue
    s1, _, _ = monthly(run, 'tsr', *LEVY)
    s0, _, _ = monthly(run, 'tsrc', *LEVY)
    a = zone(t.mean(axis=0), la) * 100
    c = zone((s1 - s0).mean(axis=0), la)
    if base is None:
        base = (a, c)
    rows.append((name, n, a, c, a - base[0], c - base[1]))
print(f'  {"run":16s} {"nyr":>4s} {"area %":>8s} {"SW CRE":>8s} {"d area":>8s} '
      f'{"d CRE":>8s}  {"area needed":>11s}')
need = o_area - base[0] if base else 0.0
for name, n, a, c, da, dc in rows:
    frac = 100 * da / need if need else 0.0
    print(f'  {name:16s} {n:4d} {a:8.2f} {c:8.2f} {da:+8.2f} {dc:+8.2f}  {frac:10.1f} %')
print(f'\n  "area needed" = fraction of the {need:+.2f} pp gap to CERES that the lever closes.')
print('  NOTE the control here is the PI arm, so its absolute area differs slightly from the')
print('  period-clean anatomy above; the DELTAS are what this table is for.')

best = max(rows[1:], key=lambda r: r[4]) if len(rows) > 1 else None
if best:
    print(f'\n  Largest area gain: {best[0]} at {best[4]:+.2f} pp of {need:+.2f} needed.')
    if abs(best[4]) < 0.1 * abs(need):
        print('  *** NO LEVER MOVES CLOUD AREA MEANINGFULLY.  D2a/D2b bought their CRE by')
        print('      making existing cloud brighter, not by making more of it.  The area')
        print('      deficit is therefore UNTOUCHED, and a round aimed at it needs a')
        print('      different knob than the INP branch -- or the INP branch pushed much')
        print('      harder than the adopted setting.')
    else:
        print('  The INP branch does move area, so pushing it further is the indicated round.')

cds.close()
