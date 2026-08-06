"""Snow STATE diagnostic: why is the soil colder when the cover is not lower?

WHY THIS EXISTS.  snow_daily_diag.py reported that O1/O2 have a POSITIVE autumn
cover anomaly (+0.14..+0.24) yet still lose 18-20 K of January soil temperature,
and concluded the autumn-seeding hypothesis was falsified.  That conclusion rested
on a reconstruction with a bug: daily() averaged sd and rsn over 40 years FIRST and
only then applied the tanh.  SCF is strongly nonlinear in depth, so
    mean_years( SCF(d) )  !=  SCF( mean_years(d) )
and the error is largest in autumn, exactly where the claim was made, because the
snow-onset date varies by weeks between years.  Everything that script said about
dCover has to be recomputed per year and only then averaged.  This script does that.

It also stops reconstructing where it does not have to.  The operator's question --
"the I series only changes snow cover, we revert that and the soil is STILL colder"
-- is answered by MEASURED state, not by a formula: if the snowpack itself (SWE,
depth, density, snow temperature) differs between the runs, then cover is not the
only channel to the soil and the whole cover-centred framing is wrong.

TWO CHANNELS, OPPOSITE SIGNS.  The source says the soil sees the snow through
    ZSNCONDH = PFRSN * (d*lambda_sn + D1*lambda_soil) / (d + D1)^2     [W m-2 K-1]
        srfsn_webal_mod.F90:239, D1 = RDAW(1) = 0.07 m,
        lambda_sn(rho) = 2.5e-6 rho^2 - 1.23e-4 rho + 0.024   (fcsurf.h:29)
and the snow layer thickness is the GRID-BOX MEAN, d = SWE/rho, NOT divided by the
snow fraction (srfsn_webal_mod.F90:222).  So
    * MORE cover  -> LARGER conductance -> soil more strongly tied to the cold pack
    * DEEPER snow -> SMALLER conductance -> insulation
Raising cover therefore COOLS the winter soil at fixed depth.  That is the reverse
of the "snow is a blanket" intuition every earlier round of this investigation used,
and it has never been tested.  The conductance below is reconstructed per day per
year from each run's own state so the two channels can be weighed against each other.

WHAT TO READ
  * cover: recomputed correctly.  Does the O1/O2 autumn surplus survive per-year
    averaging, or was it an artefact of averaging depth before applying tanh?
  * SWE / depth: if the winter pack is thinner in N2/O1/O2 the soil cools through
    conductance regardless of cover, and the cover framing is dead for a second
    reason.
  * G (conductance): the single number that should order the runs if the mechanism
    is thermal coupling.  N2 -24.6 K, O1 -19.8, O2 -17.8 in January -- G should
    rank the same way.  If it does not, the mechanism is elsewhere again.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

RUNS = [('N1 ref off', 'amip_N1_snowdiag',     None),
        ('N2 cur',     'amip_N2_snowdiag_scf', (0.016, 100.0, 1.6)),
        ('O1 m4',      'amip_O1_scf_m4',       (0.018, 170.0, 4.0)),
        ('O2 m3',      'amip_O2_scf_m3',       (0.018, 170.0, 3.0))]
YEARS = range(1876, 1916)
DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
DOY_MONTH = np.repeat(np.arange(12), DPM)

D1 = 0.07          # RDAW(1), first soil layer [m]
LAM_SOIL = 1.5     # PSURFCOND stand-in; identical across runs, so deltas are clean
RQSNCR_INV = 10.0  # as-released cover = min(1, 10*depth)


def lam_sn(rho):
    """fcsurf.h:29 -- snow thermal conductivity [W m-1 K-1]."""
    return 2.5e-6 * rho * rho - 1.23e-4 * rho + 0.024


def cover(depth, rho, p):
    """Each run's OWN cover formula: tanh branch if tuned, linear ramp if p is None."""
    if p is None:
        return np.clip(RQSNCR_INV * depth, 0.0, 1.0)
    z0, rho_new, m = p
    scale = 2.5 * z0 * (np.maximum(rho, 50.0) / rho_new) ** m
    c = np.tanh(depth / np.maximum(scale, 1e-6))
    return np.where((depth > 1e-6) & (depth * rho > 1e-6), np.clip(c, 0, 1), 0.0)


def conductance(depth, rho, frsn):
    """srfsn_webal_mod.F90:239 -- snow-to-soil conductance [W m-2 K-1]."""
    return frsn * (depth * lam_sn(rho) + D1 * LAM_SOIL) / np.maximum((depth + D1) ** 2, 1e-12)


def path(run, var, y):
    for pat in (f'atm_remapped_1d_{var}_1d_{y}-{y}.nc', f'atm_remapped_1d_{var}_{y}-{y}.nc'):
        f = f'{RT}/{run}/outdata/oifs/{pat}'
        if os.path.exists(f):
            return f
    return None


def load(run, var, y, sel):
    f = path(run, var, y)
    if f is None:
        return None
    with xr.open_dataset(f, decode_times=False) as d:
        a = d[var].values
    if a.shape[0] == 366:
        a = np.delete(a, 59, axis=0)
    return a.reshape(a.shape[0], -1)[:, sel] if a.shape[0] == 365 else None


# ---- grid and the Siberian box ------------------------------------------------
with xr.open_dataset(LSMF) as d:
    lsm = d['lsm'].isel(time_counter=0).values
    lat, lon = d['lat'].values, d['lon'].values
LA = np.broadcast_to(lat[:, None], lsm.shape)
LO = np.broadcast_to(lon[None, :], lsm.shape)
boxm = (LA >= 55) & (LA <= 75) & (LO >= 60) & (LO <= 180) & (lsm > 0.5)
sel = np.flatnonzero(boxm.ravel())
w = np.cos(np.deg2rad(LA)).ravel()[sel]
bm = lambda a: np.einsum('ij,j->i', a, w) / w.sum()      # (365, ncell) -> (365,)

print(__doc__.split('WHAT TO READ')[0])
print('=' * 92)
print(f'Siberian land 55-75N 60-180E, {sel.size} cells, per-year then averaged\n')

# ---- accumulate per year, applying every nonlinear formula BEFORE averaging ----
FIELDS = ['swe', 'depth', 'rho', 'cover', 'G', 'tsn', 'stl1', 'stl2']
S = {}
for lab, run, p in RUNS:
    acc, n = {k: np.zeros(365) for k in FIELDS}, 0
    for y in YEARS:
        sd = load(run, 'sd', y, sel)
        rsn = load(run, 'rsn', y, sel)
        if sd is None or rsn is None:
            continue
        swe = sd * 1000.0                                  # m water equiv -> kg/m2
        rho = np.maximum(rsn, 1e-6)
        depth = swe / rho
        cov = cover(depth, rho, p)
        acc['swe'] += bm(swe)
        acc['depth'] += bm(depth)
        acc['rho'] += bm(np.where(swe > 1.0, rho, np.nan) if False else rho)
        acc['cover'] += bm(cov)
        acc['G'] += bm(conductance(depth, rho, cov))
        for v in ('tsn', 'stl1', 'stl2'):
            a = load(run, v, y, sel)
            if a is not None:
                acc[v] += bm(a)
        n += 1
    if not n:
        print(f'  {lab}: no daily output'); continue
    S[lab] = {k: v / n for k, v in acc.items()}
    S[lab]['n'] = n
    print(f'  {lab:12s} {n} yr')

mon = lambda a: [a[DOY_MONTH == m].mean() for m in range(12)]
ref = S['N1 ref off']

# ---- 1. cover, recomputed correctly -------------------------------------------
print('\n\n1. SNOW COVER FRACTION (each run under its own formula, per-year averaging)\n')
print(f'  {"run":12s}' + ''.join(f'{m:>7s}' for m in MON))
for lab in S:
    print(f'  {lab:12s}' + ''.join(f'{v:7.3f}' for v in mon(S[lab]['cover'])))
print(f'\n  {"delta vs N1":12s}')
for lab in [l for l in S if l != 'N1 ref off']:
    d = np.array(mon(S[lab]['cover'])) - np.array(mon(ref['cover']))
    print(f'  {lab:12s}' + ''.join(f'{v:+7.3f}' for v in d))

# ---- 2. the pack itself: is cover really the only thing that changed? ---------
for name, unit, key in (('SNOW WATER EQUIVALENT', 'kg/m2', 'swe'),
                        ('SNOW DEPTH', 'm', 'depth'),
                        ('BULK DENSITY', 'kg/m3', 'rho')):
    print(f'\n\n2{"abc"["swe depth rho".split().index(key)]}. {name} [{unit}]\n')
    print(f'  {"run":12s}' + ''.join(f'{m:>7s}' for m in MON))
    for lab in S:
        print(f'  {lab:12s}' + ''.join(f'{v:7.2f}' for v in mon(S[lab][key])))
    print(f'  {"-- delta":12s}')
    for lab in [l for l in S if l != 'N1 ref off']:
        d = np.array(mon(S[lab][key])) - np.array(mon(ref[key]))
        print(f'  {lab:12s}' + ''.join(f'{v:+7.2f}' for v in d))

# ---- 3. the conductance -- the quantity that should order the runs -----------
print('\n\n3. SNOW-TO-SOIL CONDUCTANCE G [W m-2 K-1]   (higher = soil tied harder to cold snow)\n')
print(f'  {"run":12s}' + ''.join(f'{m:>7s}' for m in MON))
for lab in S:
    print(f'  {lab:12s}' + ''.join(f'{v:7.3f}' for v in mon(S[lab]['G'])))
print(f'  {"-- delta":12s}')
for lab in [l for l in S if l != 'N1 ref off']:
    d = np.array(mon(S[lab]['G'])) - np.array(mon(ref['G']))
    print(f'  {lab:12s}' + ''.join(f'{v:+7.3f}' for v in d))

# ---- 4. temperatures ----------------------------------------------------------
for name, key in (('SNOW TEMPERATURE tsn', 'tsn'), ('SOIL LAYER 1 stl1', 'stl1'),
                  ('SOIL LAYER 2 stl2', 'stl2')):
    print(f'\n\n4. {name} -- delta vs N1 [K]\n')
    print(f'  {"run":12s}' + ''.join(f'{m:>7s}' for m in MON))
    print(f'  {"N1 (abs)":12s}' + ''.join(f'{v:7.2f}' for v in mon(ref[key])))
    for lab in [l for l in S if l != 'N1 ref off']:
        d = np.array(mon(S[lab][key])) - np.array(mon(ref[key]))
        print(f'  {lab:12s}' + ''.join(f'{v:+7.2f}' for v in d))

# ---- 5. does G order the runs the way the soil damage does? -------------------
print('\n\n5. THE TEST -- DJF mean of each candidate against DJF soil damage\n')
djf = np.isin(DOY_MONTH, [0, 1, 11])
print(f'  {"run":12s}{"dSoil":>9s}{"dCover":>9s}{"dDepth":>9s}{"dSWE":>9s}{"dG":>9s}{"G ratio":>9s}')
for lab in [l for l in S if l != 'N1 ref off']:
    f = lambda k: S[lab][k][djf].mean() - ref[k][djf].mean()
    print(f'  {lab:12s}{f("stl1"):+9.2f}{f("cover"):+9.3f}{f("depth"):+9.3f}'
          f'{f("swe"):+9.2f}{f("G"):+9.3f}'
          f'{S[lab]["G"][djf].mean()/ref["G"][djf].mean():9.2f}')
print("""
  READING IT.  If dG ranks N2 > O1 > O2 (matching soil damage -24.6 > -19.8 > -17.8)
  then thermal coupling through the snow is the mechanism and cover was only ever a
  proxy for it.  If dCover is ~0 in DJF but dG is large, the damage is carried by
  DEPTH, not cover -- and no re-parameterisation of the cover curve can fix it.""")
