"""WHICH of the three factors in the tropical longwave CRE is deficient?

For a cloud layer the longwave cloud radiative effect is

    LW CRE  ~  f_high  x  epsilon  x  (OLR_clear - sigma T_c^4)

              area       opacity     HEIGHT the radiation escapes from

Three independent factors.  Every lever this campaign has aimed at the tropical longwave
acted on ONE of them:

    DETRPEN   sets PLUDE, the detrained condensate      -> epsilon
    RVICE     sets how fast ice falls out                -> epsilon
    RPRCON    sets how much condensate survives to detrain -> epsilon
    RLCRITSNOW / RSNOWLIN2  snow autoconversion          -> epsilon

and EPSILON SATURATES.  Once an anvil is optically thick it is already a blackbody and
more ice buys nothing.  That is a candidate physical explanation for two facts measured
in tropical_lw_scan.py: the term is extraordinarily stiff (interannual sd 0.101 W/m2,
44-yr threshold +-0.042) and yet the largest positive response among 50 arms is +0.373,
17 % of the -2.16 deficit.  If epsilon is already near 1 in the tropics, the whole
emissivity class is capped no matter which member of it we tune, and the deficit has to
be closed through AREA or HEIGHT instead -- terms nothing has ever been aimed at.

HOW THIS IS TESTED WITHOUT A NEW RUN.  Fifty arms each perturbed tropical cloud in a
different way.  Regressing d(LW CRE) on d(hcc) and d(tciw) ACROSS those arms measures the
sensitivity to area and to ice water empirically, over 50 perturbations, rather than
assuming it.  Three things fall out:

  * the tciw slope, and whether it FLATTENS at high tciw -- the saturation test
  * the hcc slope, hence how much area a +2.16 closure would actually need
  * the RESIDUAL, which is what neither area nor ice water explains.  Since those are two
    of the three factors, the residual IS the cloud-top-height signal.

WHY NOT JUST COMPARE TO CERES.  Because of the trap round 23 paid for: CERES cldarea is a
MODIS *mask* that counts optically thin cloud, model hcc is a *radiative* cover, and the
two are not the same quantity; and ERA5 is IFS, so it is not an independent reference for
an IFS-derived model.  This regression is model-internal.  The only observational number
it uses is the LW CRE flux deficit itself, which is a radiative flux and IS directly
comparable.  CERES cldtau and cldtemp are printed alongside as a cross-check, with that
definitional caveat attached rather than hidden.

PERIOD.  All arms PI-epoch over the campaign-standard 1872-1915, deltas against
amip_pi_base, so every point is epoch-clean.  Run list comes from runs.py.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import sys
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runs import RUNS, RT, Y0, Y1

ACC = 3600.0
CTRL = 'amip_pi_base'
TROP = (-30.0, 30.0)
DEFICIT = 2.16                 # W/m2 of tropical LW CRE still to be found
K_ICE = 0.06                   # m2/g, LW mass absorption coefficient for cloud ice
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
VARS = ('ttr', 'ttrc', 'hcc', 'tciw', 'tcc', 'tclw')

print(__doc__)
print('=' * 104)


def load(run):
    """Tropics-mean of each field, averaged over Y0..Y1.  None if the run is short."""
    D = f'{RT}/{run}/outdata/oifs'
    acc = {v: [] for v in VARS}
    lat = None
    for y in range(Y0, Y1 + 1):
        paths = {v: f'{D}/atm_remapped_1m_{v}_1m_{y}-{y}.nc' for v in VARS}
        if not all(os.path.exists(p) for p in paths.values()):
            return None
        for v, p in paths.items():
            with xr.open_dataset(p, decode_times=False) as ds:
                a = ds[v].values
                if lat is None:
                    lat = ds['lat'].values
            acc[v].append(a / (ACC if v in ('ttr', 'ttrc') else 1.0))
    sel = (lat >= TROP[0]) & (lat < TROP[1])
    w = np.cos(np.deg2rad(lat[sel]))

    def bm(a4):
        m = np.mean(a4, axis=0).mean(axis=0)        # years, then months
        return float(np.average(m[sel, :].mean(axis=1), weights=w))
    out = {v: bm(np.array(acc[v])) for v in VARS}
    out['lwcre'] = out['ttr'] - out['ttrc']
    return out


ctl = load(CTRL)
if ctl is None:
    raise SystemExit(f'{CTRL} does not cover {Y0}-{Y1}')

rows = []
for label, r in RUNS:
    if r == CTRL:
        continue
    d = load(r)
    if d is None:
        continue
    rows.append((label, d))

print(f'control {CTRL}, {Y0}-{Y1}; {len(rows)} arms scored\n')
print('1. THE CONTROL STATE -- is the tropical anvil already a blackbody?')
print('-' * 104)
iwp_box = ctl['tciw'] * 1000.0                       # grid-box mean, g/m2
iwp_in = iwp_box / max(ctl['hcc'], 1e-6)             # in-cloud, dividing by high-cloud area
eps_box = 1.0 - np.exp(-K_ICE * iwp_box)
eps_in = 1.0 - np.exp(-K_ICE * iwp_in)
print(f'  tropical high cloud cover  hcc  = {ctl["hcc"]:.4f}')
print(f'  total cloud cover          tcc  = {ctl["tcc"]:.4f}')
print(f'  column ice water          tciw  = {ctl["tciw"]*1000:.2f} g/m2  (grid-box mean)')
print(f'  column liquid water       tclw  = {ctl["tclw"]*1000:.2f} g/m2')
print(f'  LW CRE                          = {ctl["lwcre"]:.3f} W/m2   (deficit {DEFICIT:+.2f})')
print(f'\n  emissivity from grid-box IWP, eps = 1-exp(-{K_ICE}*IWP): {eps_box:.3f}')
print(f'  emissivity from IN-CLOUD IWP (IWP/hcc = {iwp_in:.1f} g/m2): {eps_in:.3f}')
print('  (in-cloud is the physically relevant one; grid-box mean understates it)')
if eps_in > 0.9:
    print('  => SATURATED.  The anvil is already effectively a blackbody, so every')
    print('     emissivity lever (DETRPEN, RVICE, RPRCON, RLCRITSNOW) is capped.')
elif eps_in > 0.7:
    print('  => PARTLY saturated: real but strongly diminishing headroom in emissivity.')
else:
    print('  => NOT saturated: emissivity levers still have room, so the 17 % ceiling')
    print('     measured in the scan needs a different explanation.')

# ------------------------------------------------------------------ 2. the regression
print('\n2. REGRESSION ACROSS THE ARCHIVE: d(LW CRE) = a*d(hcc) + b*d(tciw)')
print('-' * 104)
dy = np.array([d['lwcre'] - ctl['lwcre'] for _, d in rows])
dA = np.array([d['hcc'] - ctl['hcc'] for _, d in rows])
dI = np.array([(d['tciw'] - ctl['tciw']) * 1000.0 for _, d in rows])
X = np.column_stack([dA, dI])
coef, *_ = np.linalg.lstsq(X, dy, rcond=None)
pred = X @ coef
resid = dy - pred
ss = 1.0 - np.sum(resid**2) / np.sum(dy**2)
print(f'  d(LW CRE) = {coef[0]:+.3f} * d(hcc)  {coef[1]:+.4f} * d(tciw[g/m2])'
      f'      R2 = {ss:.3f}  (n={len(dy)})')
print(f'  area  slope: {coef[0]:+.3f} W/m2 per unit hcc fraction '
      f'(= {coef[0]/100:+.4f} per pp)')
print(f'  ice   slope: {coef[1]:+.4f} W/m2 per g/m2 of column ice')
print(f'\n  to buy {DEFICIT:+.2f} W/m2 through AREA alone: '
      f'd(hcc) = {DEFICIT/coef[0]*100:+.2f} pp'
      if abs(coef[0]) > 1e-9 else '')
print(f'  to buy {DEFICIT:+.2f} W/m2 through ICE  alone: '
      f'd(tciw) = {DEFICIT/coef[1]:+.1f} g/m2 '
      f'({DEFICIT/coef[1]/max(iwp_box,1e-9)*100:+.0f} % of the control column)'
      if abs(coef[1]) > 1e-9 else '')

print(f'\n  variance explained by area+ice: {100*ss:.1f} %'
      f'   -> residual (the HEIGHT term) carries {100*(1-ss):.1f} %')

print('\n3. THE ARMS THE TWO-FACTOR MODEL CANNOT EXPLAIN  (largest residuals)')
print('-' * 104)
print(f'  {"run":18s} {"d LWCRE":>9s} {"predicted":>10s} {"residual":>9s} '
      f'{"d hcc [pp]":>11s} {"d tciw":>9s}')
order = np.argsort(-np.abs(resid))
for k in order[:12]:
    print(f'  {rows[k][0]:18s} {dy[k]:9.3f} {pred[k]:10.3f} {resid[k]:9.3f} '
          f'{dA[k]*100:11.3f} {dI[k]:9.3f}')
print('\n  A large positive residual = more LW CRE than its area and ice water explain,')
print('  i.e. the cloud moved UP.  Those arms are where the height term is visible.')

# ------------------------------------------------------------------ 3. CERES cross-check
print('\n4. CERES CROSS-CHECK (definitional caveat applies -- MODIS mask, not radiative cover)')
print('-' * 104)
try:
    with xr.open_dataset(CERESF) as c:
        clat = c['lat'].values
        s = (clat >= TROP[0]) & (clat < TROP[1])
        w = np.cos(np.deg2rad(clat[s]))

        def cz(v):
            a = c[v].values.mean(axis=0)
            return float(np.average(a[s, :].mean(axis=1), weights=w))
        for v, unit in (('cldarea_total_daynight_clim', '%'),
                        ('cldpress_total_daynight_clim', 'hPa'),
                        ('cldtemp_total_daynight_clim', 'K'),
                        ('cldtau_total_day_clim', '-')):
            if v in c.data_vars:
                print(f'  CERES {v:38s} {cz(v):9.2f} {unit}')
    print('\n  cldpress/cldtemp are the observational handle on the HEIGHT term; the model')
    print('  writes no cloud-top pressure, so this is the only external check available,')
    print('  and it is an all-cloud effective value, not a high-cloud one.')
except FileNotFoundError:
    print('  CERES file not found -- skipped')
