"""Where does the March-April melt energy go? Close the surface energy budget.

snow_budget.py established three things: spring snowfall is right, the loss term is
mistimed (Mar -13.7, Apr -14.0, then May +29.3 mm/month vs ERA5), and it is NOT an
energy-supply problem -- net surface SW is POSITIVE against CERES in exactly the
months that fail to melt. So the energy is present and the snowpack is not turning
it into melt. This script asks where it goes instead.

Surface energy balance over land, all terms in W/m2, IFS sign convention (fluxes
POSITIVE DOWNWARD, i.e. into the surface):

    ssr + str + sshf + slhf  =  G + M

    ssr   net shortwave        str   net longwave
    sshf  sensible heat        slhf  latent heat (includes sublimation)
    G     conduction into snow/soil, plus the energy consumed warming the pack
    M     energy consumed melting  =  L_f * melt_rate,  L_f = 3.34e5 J/kg

Three hypotheses this separates:

  1. COLD CONTENT -- the pack is too cold, so absorbed energy warms it toward 0 C
     instead of melting it. Signature: model `tsn` colder than ERA5 in Feb-Apr, and
     a larger residual G.
  2. TURBULENT LOSS -- the absorbed energy leaves as sensible/latent flux to the
     atmosphere rather than melting. Signature: model sshf/slhf more negative
     (more upward) than ERA5.
  3. LONGWAVE LOSS -- the surface radiates the gain away. Signature: model `str`
     more negative than ERA5.

ERA5 provides every term (146 sshf, 147 slhf, 176 ssr, 177 str, 238 tsn), so this is
a like-for-like comparison rather than a model-only diagnosis.

UNITS, asserted not assumed (a previous pass was wrong by 10^3):
  * model accumulated fluxes  : / ACC (3600 s) -> W/m2
  * ERA5 fc monthly means     : accumulated per DAY -> / 86400 -> W/m2
  * melt energy               : mm/month w.e. -> kg/m2/month -> * L_f / seconds -> W/m2
A printed closure check reports the residual; if it is comparable to the terms
themselves the decomposition is not to be trusted.
"""
import numpy as np, xarray as xr, os, sys, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF, Y0, Y1

BOX = ((55, 75), (60, 180))
ACC = 3600.0
LF = 3.34e5                                   # latent heat of fusion, J/kg
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
DPM = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
SPM = DPM * 86400.0
E5 = '/work/ab0246/a270092/obs/era5/snow'
PD = list(range(1990, 2015))
PI = list(range(Y0, Y1 + 1))

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def bm(a, lat, lon):
    ys = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
    l180 = ((lon + 180) % 360) - 180
    xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub = a[ii]
    m = np.isfinite(sub) & (lsm[ii] > 0.5)
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape)
    return np.average(sub[m], weights=w[m]) if m.any() else np.nan


def clim(run, var, yrs):
    acc, n, lat, lon = None, 0, None, None
    for y in yrs:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        d = xr.open_dataset(f)
        a = d[var].values
        lat, lon = d[var].lat.values, d[var].lon.values
        d.close()
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, lat, lon) if n else (None, None, None)


# ------------------------------------------------------------------ ERA5
_d = xr.open_dataset(os.path.join(E5, 'era5_lsm.nc'))
lsm5 = np.squeeze(_d['var172'].values); _d.close()


def e5(pcode, vn):
    p = os.path.join(E5, f'era5_{pcode}_clim_1990-2014.nc')
    if not os.path.exists(p):
        return None, None, None
    d = xr.open_dataset(p)
    a = d[vn].values
    la, lo = d['lat'].values, d['lon'].values
    d.close()
    return a, la, lo


def bm5(a, la, lo):
    ys = (la >= BOX[0][0]) & (la <= BOX[0][1])
    l180 = ((lo + 180) % 360) - 180
    xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub = a[ii]
    m = np.isfinite(sub) & (lsm5[ii] > 0.5)
    w = np.broadcast_to(np.cos(np.deg2rad(la[ys]))[:, None], sub.shape)
    return np.average(sub[m], weights=w[m]) if m.any() else np.nan


E = {}
for code, vn, key in (('176', 'var176', 'ssr'), ('177', 'var177', 'str'),
                      ('146', 'var146', 'sshf'), ('147', 'var147', 'slhf'),
                      ('045', 'var45', 'melt'), ('144', 'var144', 'sf'),
                      ('141', 'var141', 'sd'), ('238', 'var238', 'tsn')):
    a, la, lo = e5(code, vn)
    if a is None:
        sys.exit(f'ERA5 {code} missing -- rerun albedo_decompose_prep.sh')
    v = np.array([bm5(a[m], la, lo) for m in range(12)])
    if key in ('ssr', 'str', 'sshf', 'slhf'):
        v = v / 86400.0                       # J/m2/day -> W/m2
    elif key in ('melt', 'sf'):
        v = v * 1000.0 * DPM                  # m/day -> mm/month
    elif key == 'sd':
        v = v * 1000.0                        # m -> mm
    E[key] = v
E['meltW'] = E['melt'] / 1000.0 * 1000.0 * LF / SPM      # mm/month -> W/m2
E['res'] = E['ssr'] + E['str'] + E['sshf'] + E['slhf'] - E['meltW']


# ------------------------------------------------------------------ model
def model(run, yrs):
    M = {}
    for v in ('ssr', 'str', 'sshf', 'slhf', 'sf', 'sd', 'tsn'):
        a, lat, lon = clim(run, v, yrs)
        if a is None:
            return None
        x = np.array([bm(a[m], lat, lon) for m in range(12)])
        if v in ('ssr', 'str', 'sshf', 'slhf'):
            x = x / ACC
        elif v == 'sf':
            x = x / ACC * SPM * 1000.0
        elif v == 'sd':
            x = x * 1000.0
        M[v] = x
    dsd = np.roll(M['sd'], -1) - M['sd']
    M['loss'] = M['sf'] - dsd                              # melt + sublimation, mm/month
    M['lossW'] = M['loss'] * LF / SPM                      # as W/m2 if it were all melt
    M['res'] = M['ssr'] + M['str'] + M['sshf'] + M['slhf'] - M['lossW']
    return M


runs = [('presentday', 'amip_presentday', PD), ('control', 'amip_pi_base', PI),
        ('G4 tundra', 'amip_G4_tundra', PI)]
R = {}
for lab, r, yrs in runs:
    m = model(r, yrs)
    if m is None:
        print(f'  !! {lab} missing'); continue
    R[lab] = m

MM = [1, 2, 3, 4]                                          # Feb Mar Apr May

print('=' * 96)
print('SURFACE ENERGY BUDGET over Siberian land [W/m2], IFS sign convention (positive DOWN)')
print('period-clean: model amip_presentday (1990-2014) vs ERA5 (1990-2014)')
print('=' * 96)
ref = R.get('presentday')
print(f'\n  {"":6s} {"ssr":>17s} {"str":>17s} {"sshf":>17s} {"slhf":>17s}')
print(f'  {"":6s} ' + ' '.join(f'{"model":>8s} {"ERA5":>8s}' for _ in range(4)))
for m in MM:
    row = ''
    for k in ('ssr', 'str', 'sshf', 'slhf'):
        row += f' {ref[k][m]:8.1f} {E[k][m]:8.1f}'
    print(f'  {MON[m]:6s}{row}')

print(f'\n  {"":6s} {"Rnet+turb":>19s} {"melt energy":>19s} {"residual (G)":>19s}')
print(f'  {"":6s} ' + ' '.join(f'{"model":>9s} {"ERA5":>9s}' for _ in range(3)))
for m in MM:
    am = ref['ssr'][m] + ref['str'][m] + ref['sshf'][m] + ref['slhf'][m]
    ae = E['ssr'][m] + E['str'][m] + E['sshf'][m] + E['slhf'][m]
    print(f'  {MON[m]:6s} {am:9.1f} {ae:9.1f} {ref["lossW"][m]:9.1f} {E["meltW"][m]:9.1f} '
          f'{ref["res"][m]:9.1f} {E["res"][m]:9.1f}')

print('\n  (model "melt energy" uses the COMBINED loss, so it is an upper bound --')
print('   part of it is sublimation, which is already counted in slhf. ERA5 uses its')
print('   reported melt field. The residual absorbs conduction + pack warming.)')

print('\n' + '=' * 96)
print('HYPOTHESIS 1 -- COLD CONTENT: is the model snowpack too cold?  `tsn` [K]')
print('=' * 96)
print(f'  {"":6s} {"ERA5":>9s} ' + ' '.join(f'{l:>14s}' for l in R))
for m in range(12):
    row = ' '.join(f'{R[l]["tsn"][m]-E["tsn"][m]:>+14.2f}' for l in R)
    mk = '  <<<' if m in MM else ''
    print(f'  {MON[m]:6s} {E["tsn"][m]:9.2f} {row}{mk}')
print('  columns are model MINUS ERA5 snow temperature')

print('\n' + '=' * 96)
print('HYPOTHESIS 2/3 -- where the energy goes: model MINUS ERA5 [W/m2]')
print('=' * 96)
print(f'  {"":6s} {"ssr":>8s} {"str":>8s} {"sshf":>8s} {"slhf":>8s} {"avail":>8s} {"meltE":>8s} {"resid":>8s}')
for m in MM:
    am = ref['ssr'][m] + ref['str'][m] + ref['sshf'][m] + ref['slhf'][m]
    ae = E['ssr'][m] + E['str'][m] + E['sshf'][m] + E['slhf'][m]
    print(f'  {MON[m]:6s} {ref["ssr"][m]-E["ssr"][m]:8.1f} {ref["str"][m]-E["str"][m]:8.1f} '
          f'{ref["sshf"][m]-E["sshf"][m]:8.1f} {ref["slhf"][m]-E["slhf"][m]:8.1f} '
          f'{am-ae:8.1f} {ref["lossW"][m]-E["meltW"][m]:8.1f} {ref["res"][m]-E["res"][m]:8.1f}')
print('\n  A negative sshf/slhf difference = model loses MORE to the atmosphere.')
print('  A negative str difference        = model radiates MORE away.')
print('  A positive residual difference   = model puts MORE into warming/conducting.')


# ---------------------------------------------------------------------------
# ⚠ CAVEAT ON THIS SCRIPT'S CENTRAL COMPARISON (2026-08-03).
#
# The turbulent-flux terms here (sshf, slhf) are compared against ERA5, but ERA5's
# surface fluxes are PURE HTESSEL OUTPUT -- the same land-surface family as ours.
# No observation constrains them. The apparent result (model vents ~3-5 W/m2 more
# upward sensible heat than ERA5 in Feb-Apr) therefore compares our scheme to a
# sibling of itself and carries NO observational authority. Treat it as a
# hypothesis to be tested, never as a finding.
#
# What IS observationally supported, from CERES (a radiative retrieval) and the
# Rutgers/IMS satellite snow record (see snowcover_vs_satellite.py):
#   * surface albedo: model +0.046 too high in June, +0.037 in May, and -0.049 /
#     -0.040 too LOW in Mar/Apr.  ERA5 agrees with CERES to +-0.02 in all sunlit
#     months, which is what validates ERA5's albedo for this purpose.
#   * snow cover: model June fraction 0.380 vs satellite 0.203, and melt-out
#     (50 % crossing) 24 May vs satellite 11 May -- about 13 days late.
#
# The cold-content hypothesis is separately REJECTED on ERA5 `tsn` (model within a
# few tenths of a K), but that too is a model-vs-model comparison.
#
# Only flux towers can settle the turbulent-flux question. Siberian candidates are
# point sites (Spasskaya Pad, ZOTTO), not a gridded product.
# ---------------------------------------------------------------------------
