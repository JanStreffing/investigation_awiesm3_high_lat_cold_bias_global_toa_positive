"""Does anything in the archive already move the TROPICAL LONGWAVE cloud term?

WHY THIS, AND WHY NOW.  The four-way attribution (round 22) found the two remaining
regional failures act through DIFFERENT channels:

    Southern Ocean 65-45S   fails through SW cloud   +7.45 W/m2  (too little reflection)
    tropics 30S-30N         fails through LW cloud   -2.16 W/m2  (too little trapping)

The Southern Ocean side now has a lever that works -- DMS-derived CCN, which buys 67 %
of the DJF deficit -- but it is blocked because it also brightens the tropics, and the
tropics cannot afford to lose any more absorbed energy.  So the tropical LW deficit is
not a separate item on the list any more: it GATES the Southern Ocean fix.  Close it and
the DMS lever becomes affordable.

There is also a global-budget reason to want the pair rather than either alone.  They
push opposite ways: DMS CCN removes energy globally (-1.72 W/m2 at S=166), which helps
the positive imbalance, while adding tropical LW trapping puts energy back.  Run
together they fix both regional biases and partly cancel in the global mean.

WHAT WOULD MAKE A LEVER USABLE HERE.  Not "acts through the longwave" -- that is too
loose.  RVICE acts through the longwave and is exactly the counter-example: it is an ice
fall speed, it applies wherever there is cloud ice, and it therefore TRADES Southern
Ocean high cloud against tropical anvils.  It was raised 0.13 -> 0.16 to suppress SO hcc
(project_management #169/#170), at a cost to the tropics.  A usable lever has to be keyed
on DEEP CONVECTION, which barely exists over the Southern Ocean, so that the selectivity
comes from where convection happens rather than from a latitude band.

WHY A SCAN BEFORE A RUN.  Roughly 45 % of this campaign's retracted claims were
answerable from data already on disk.  Fifty-odd runs have been scored on Southern Ocean
SW and on global RMSE; none has been scored on the tropical LW term.  If something in the
archive already moves it, that is a 48-year run we do not have to spend -- and if
nothing does, that null is itself the finding, and it tells us the lever has to come from
outside the parameters tried so far.

METHOD.  All arms are PI-epoch, evaluated over the campaign-standard 1872-1915 (44 yr),
against the shared control amip_pi_base over the same years, so every delta is
epoch-clean by construction.  Detection thresholds are computed FIRST, from the control's
own interannual scatter, as 1.96 * sd * sqrt(2/44) for a 44-year pair -- the same
construction that gives the known +-1.97 for a 1-year Southern Ocean SW pair.

Signs: LW CRE = ttr - ttrc is POSITIVE (cloud reduces outgoing longwave).  We want the
tropical value to go UP.  SW CRE = tsr - tsrc is NEGATIVE; we want the Southern Ocean
value to go DOWN (more negative = more reflection).  Net TOA = tsr + ttr.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

ACC = 3600.0
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
CTRL = 'amip_pi_base'
Y0, Y1 = 1872, 1915                      # campaign standard, 44 years, PI epoch
TROP, SO = (-30.0, 30.0), (-65.0, -45.0)
TARGET_TROP_LW = +2.16                   # what closing the tropical deficit requires
# Withdrawn in round 23 -- kept in the table but marked, never silently dropped.
WITHDRAWN = {'amip_S1_clddiff3e6', 'amip_S2_ovlliqice035', 'amip_S3_clcritsea6e4'}

print(__doc__)
print('=' * 108)


def load(run):
    """Per-year band means for the four fields; None if the run lacks the full period."""
    D = f'{RT}/{run}/outdata/oifs'
    years = []
    for y in range(Y0, Y1 + 1):
        f = {v: f'{D}/atm_remapped_1m_{v}_1m_{y}-{y}.nc' for v in
             ('tsr', 'tsrc', 'ttr', 'ttrc')}
        if not all(os.path.exists(p) for p in f.values()):
            return None, None
        d = {}
        for v, p in f.items():
            with xr.open_dataset(p, decode_times=False) as ds:
                d[v] = ds[v].values / ACC
                lat = ds['lat'].values
        if d['tsr'].shape[0] != 12:
            return None, None
        years.append((d, lat))
    return years, years[0][1]


def band(a3d, lat, lo, hi):
    sel = (lat >= lo) & (lat < hi)
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(a3d.mean(axis=0)[sel, :].mean(axis=1), weights=w))


def metrics(years, lat):
    """(nyear, 5) array: trop LW CRE, trop SW CRE, SO SW CRE, SO LW CRE, global net TOA."""
    out = []
    for d, _ in years:
        swcre, lwcre = d['tsr'] - d['tsrc'], d['ttr'] - d['ttrc']
        net = d['tsr'] + d['ttr']
        out.append([band(lwcre, lat, *TROP), band(swcre, lat, *TROP),
                    band(swcre, lat, *SO), band(lwcre, lat, *SO),
                    band(net, lat, -90, 90)])
    return np.array(out)


NAMES = ['trop LW CRE', 'trop SW CRE', 'SO SW CRE', 'SO LW CRE', 'global net TOA']

cy, lat = load(CTRL)
if cy is None:
    raise SystemExit(f'control {CTRL} does not cover {Y0}-{Y1}')
cm = metrics(cy, lat)
n = len(cm)
thr = 1.96 * cm.std(axis=0, ddof=1) * np.sqrt(2.0 / n)

print(f'control {CTRL}, {n} years {Y0}-{Y1}\n')
print('0. DETECTION THRESHOLDS, computed before any comparison')
print('-' * 108)
print(f'  {"metric":18s} {"control mean":>13s} {"sd(interannual)":>17s} '
      f'{"44-yr pair 95%":>16s}')
for i, nm in enumerate(NAMES):
    print(f'  {nm:18s} {cm[:, i].mean():13.3f} {cm[:, i].std(ddof=1):17.3f} '
          f'{thr[i]:16.3f}')
print(f'\n  the tropical LW deficit to be closed is {TARGET_TROP_LW:+.2f} W/m2, i.e. '
      f'{TARGET_TROP_LW / thr[0]:.0f}x the detection threshold')

runs = sorted(r for r in os.listdir(RT)
              if r.startswith('amip_') and r != CTRL and
              os.path.isdir(f'{RT}/{r}/outdata/oifs'))
rows, skipped = [], []
for r in runs:
    ry, rlat = load(r)
    if ry is None:
        skipped.append(r)
        continue
    d = metrics(ry, rlat).mean(axis=0) - cm.mean(axis=0)
    rows.append((r, d))

print(f'\n1. EVERY 44-YEAR ARM, RANKED BY TROPICAL LW CRE RESPONSE   '
      f'({len(rows)} scored, {len(skipped)} lack {Y0}-{Y1})')
print('-' * 108)
print(f'  {"run":32s} ' + ' '.join(f'{nm:>15s}' for nm in NAMES) + '  flags')
rows.sort(key=lambda t: -t[1][0])
hits = []
for r, d in rows:
    f = []
    if abs(d[0]) > thr[0]:
        f.append('trop LW MOVES')
        if d[0] > 0:
            hits.append((r, d))
    if abs(d[2]) > thr[2]:
        f.append('SO SW side-effect')
    if r in WITHDRAWN:
        f.append('WITHDRAWN r23')
    print(f'  {r:32s} ' + ' '.join(f'{v:15.3f}' for v in d) + '  ' + '; '.join(f))
if skipped:
    print(f'\n  not scored (no {Y0}-{Y1}): ' + ', '.join(skipped))

print('\n2. THE ANSWER')
print('-' * 108)
if not hits:
    print('  NOTHING in the archive raises the tropical LW cloud term above its own\n'
          '  detection threshold.  The lever has to come from outside the parameters\n'
          '  tried so far, and by the argument above it should be keyed on deep\n'
          '  convection so that it cannot reach the Southern Ocean.')
else:
    print(f'  {len(hits)} arm(s) raise the tropical LW term.  Judge each on whether the\n'
          '  Southern Ocean pays for it -- a lever that moves both is RVICE again:\n')
    for r, d in hits:
        frac = 100 * d[0] / TARGET_TROP_LW
        clean = ('SO untouched' if abs(d[2]) <= thr[2]
                 else f'but SO SW moves {d[2]:+.2f} (threshold {thr[2]:.2f})')
        print(f'    {r:32s} trop LW {d[0]:+.3f} ({frac:.0f} % of the deficit), {clean}')
