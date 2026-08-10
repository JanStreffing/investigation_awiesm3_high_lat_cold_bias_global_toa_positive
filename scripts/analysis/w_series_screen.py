"""The v3.4.2 DARS2 cloud stack, taken apart and put back together, at one year.

WHAT THIS IS.  awiesm3-v3.4.2-levante-TCO95L91-DARS2_{1d,2y}.yaml carry a jointly-tuned
eight-parameter cloud stack for our exact resolution:

    namcumf: RPRCON 1.48E-3   ENTRORG 2.07E-3   RMFDEPS 0.48   DETRPEN 1.32E-4
             ENTRDD 1.08E-4
    namcldp: RVICE 0.18       RLCRITSNOW 1.46E-5   RSNOWLIN2 0.04

No fort.4 on disk sets any of them, so as far as this machine knows the stack was written
and never run.  The campaign has meanwhile been testing one knob at a time away from a
defaults base, without ever seeing it -- preflight.py was scanning a single runtime tree
and a single runscript directory.  Fixed 2026-08-10, after the same blind spot hid the
RPRCON history (#87/#95) earlier the same day.

WHY THE COMBINATION IS THE POINT.  These were optimised together, and superposition has
been wrong IN SIGN twice in this campaign (AB, ABB8).  W9 minus the sum of W1..W8 measures
the interaction directly rather than assuming it away.  A large interaction term means the
individual arms cannot be read as a menu.

WHAT ONE YEAR LICENSES.  Only the radiative metrics, and only those above their own
thresholds, computed from the control's interannual scatter over 1872-1915:

    tropical LW CRE  +-0.279     SO SW CRE  +-1.974
    tropical SW CRE  +-0.990     SO LW CRE  +-0.584
    tropical net TOA +-1.656     global net TOA +-1.333

Anything below its threshold is a BOUND, not a result.  Nothing here licenses a
temperature claim -- the seasonal T2m thresholds are 44-year numbers (DJF +-0.588,
JJA +-0.242) and a single year cannot approach them.

TARGETS, for reading the signs.  Tropical LW CRE is -2.16 short (we want it UP).  SO SW
CRE is +7.85 short of CERES (we want it DOWN, more negative).  Global net TOA is +0.64
against a PI target near zero (we want it DOWN).  The tropics sit only -0.67 from CERES
period-clean, so a tropical net TOA move beyond ~0.5 is a global knob in disguise.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
ACC, YEAR, CTRL = 3600.0, 1870, 'amip_pi_base'
ARMS = [('W1 RPRCON 1.48e-3', 'amip_W1_rprcon148'),
        ('W2 ENTRORG 2.07e-3', 'amip_W2_entrorg207'),
        ('W3 RMFDEPS 0.48', 'amip_W3_rmfdeps048'),
        ('W4 DETRPEN 1.32e-4', 'amip_W4_detrpen132'),
        ('W5 ENTRDD 1.08e-4', 'amip_W5_entrdd108'),
        ('W6 RVICE 0.18', 'amip_W6_rvice018'),
        ('W7 RLCRITSNOW 1.46e-5', 'amip_W7_lcritsnow146'),
        ('W8 RSNOWLIN2 0.04', 'amip_W8_rsnowlin204')]
COMBO = ('W9 DARS2 stack', 'amip_W9_dars2stack')
# measured from amip_pi_base 1872-1915, 1.96*sd*sqrt(2) for a pair of single years
THR = {'trop LW CRE': 0.279, 'trop SW CRE': 0.990, 'trop net TOA': 1.656,
       'SO SW CRE': 1.974, 'SO LW CRE': 0.584, 'global net TOA': 1.333}
KEYS = list(THR)

print(__doc__)
print('=' * 112)


def metrics(run):
    D = f'{RT}/{run}/outdata/oifs'
    d = {}
    for v in ('tsr', 'tsrc', 'ttr', 'ttrc'):
        f = f'{D}/atm_remapped_1m_{v}_1m_{YEAR}-{YEAR}.nc'
        if not os.path.exists(f):
            return None
        with xr.open_dataset(f, decode_times=False) as ds:
            d[v] = ds[v].values / ACC
            lat = ds['lat'].values
    sw, lw, net = d['tsr'] - d['tsrc'], d['ttr'] - d['ttrc'], d['tsr'] + d['ttr']

    def b(a, lo, hi):
        s = (lat >= lo) & (lat < hi)
        w = np.cos(np.deg2rad(lat[s]))
        return float(np.average(a.mean(axis=0)[s, :].mean(axis=1), weights=w))
    return {'trop LW CRE': b(lw, -30, 30), 'trop SW CRE': b(sw, -30, 30),
            'trop net TOA': b(net, -30, 30), 'SO SW CRE': b(sw, -65, -45),
            'SO LW CRE': b(lw, -65, -45), 'global net TOA': b(net, -90, 90)}


c = metrics(CTRL)
if c is None:
    raise SystemExit(f'control {CTRL} year {YEAR} missing')
print(f'control {CTRL}, year {YEAR}\n')
print(f'  {"":24s} ' + ' '.join(f'{k:>15s}' for k in KEYS))
print(f'  {"CONTROL (absolute)":24s} ' + ' '.join(f'{c[k]:15.3f}' for k in KEYS))
print(f'  {"1-yr pair threshold":24s} ' + ' '.join(f'{THR[k]:15.3f}' for k in KEYS))
print('-' * 112)

deltas, missing = {}, []
for lab, r in ARMS + [COMBO]:
    m = metrics(r)
    if m is None:
        missing.append(lab)
        continue
    d = {k: m[k] - c[k] for k in KEYS}
    deltas[lab] = d
    flag = ''.join('*' if abs(d[k]) > THR[k] else ' ' for k in KEYS)
    print(f'  {lab:24s} ' + ' '.join(f'{d[k]:15.3f}' for k in KEYS) + f'   {flag}')
print("\n  '*' marks a value above its own 1-year threshold; unmarked entries are bounds.")
if missing:
    print(f'  not yet on disk: {", ".join(missing)}')

# ------------------------------------------------------------------ interaction
if COMBO[0] in deltas and all(l in deltas for l, _ in ARMS):
    print('\nINTERACTION: does the stack equal the sum of its parts?')
    print('-' * 112)
    print(f'  {"":24s} ' + ' '.join(f'{k:>15s}' for k in KEYS))
    s = {k: sum(deltas[l][k] for l, _ in ARMS) for k in KEYS}
    combo = deltas[COMBO[0]]
    inter = {k: combo[k] - s[k] for k in KEYS}
    print(f'  {"sum of W1..W8":24s} ' + ' '.join(f'{s[k]:15.3f}' for k in KEYS))
    print(f'  {"W9 measured":24s} ' + ' '.join(f'{combo[k]:15.3f}' for k in KEYS))
    print(f'  {"INTERACTION (W9 - sum)":24s} ' + ' '.join(f'{inter[k]:15.3f}' for k in KEYS)
          + '   ' + ''.join('*' if abs(inter[k]) > THR[k] else ' ' for k in KEYS))
    big = [k for k in KEYS if abs(inter[k]) > THR[k]]
    if big:
        print(f'\n  NON-ADDITIVE on: {", ".join(big)}.  The arms are not a menu -- the stack'
              '\n  has to be judged, and adopted or rejected, as one object.')
    else:
        print('\n  Additive to within the thresholds: the arms CAN be read individually,'
              '\n  and a subset can be assembled without re-running the combination.')

# ------------------------------------------------------------------ verdict
if COMBO[0] in deltas:
    d = deltas[COMBO[0]]
    print('\nTHE STACK AGAINST THE CAMPAIGN TARGETS')
    print('-' * 112)
    print(f'  SO SW CRE      {d["SO SW CRE"]:+7.3f}  (want negative; gap to CERES is 7.85)'
          f'   -> {100*max(-d["SO SW CRE"],0)/7.85:.0f} % of the gap')
    print(f'  tropical LW    {d["trop LW CRE"]:+7.3f}  (want positive; deficit 2.16)'
          f'         -> {100*max(d["trop LW CRE"],0)/2.16:.0f} % of the deficit')
    print(f'  global net TOA {d["global net TOA"]:+7.3f}  (want negative; control +0.64)')
    print(f'  tropical net   {d["trop net TOA"]:+7.3f}  (guardrail ~0.5; tropics only '
          f'-0.67 from CERES)')
