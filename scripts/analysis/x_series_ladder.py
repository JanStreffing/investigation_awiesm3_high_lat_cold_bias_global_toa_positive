"""The X ladder: everything adopted, plus the two new levers, one step at a time.

    T3   adopted stack, DMS field read, ECE_DMS_CCN_SENS = 0     <- control
    X1   + RSNOWLIN2 0.04                                        (W8 on the real base)
    X2   + ECE_DMS_CCN_SENS 166                                  (U1 on top of W8)
    X3   + RCL_INPPMIN 50000                                     (S4 on top of that)

The adopted stack is G4 (RVRSMIN 1000/1000/225), K1 (RVVEGALB types 1/10/11 and SOILALB
x0.95), P3 fitted snow depletion at SWEMIN=15, D2b (RCL_INPSEA 0.2, RCL_INPPMIN 70000),
RVICE 0.16, and the DMS-Rev3 ICMCL.

WHY A LADDER.  W9 combined the eight DARS2 parameters and came out significantly
NON-ADDITIVE -- the Southern Ocean longwave interaction was +2.415 W/m2, with the sum of
the parts predicting -1.561 and the combination measuring +0.854.  Superposition has also
been wrong in SIGN twice before in this campaign (AB, ABB8).  So the interesting number
here is not the total; it is whether each lever still does on this base what it did alone.

THE COMPARISON THAT MATTERS.  Each rung is scored against its own control AND against the
same lever measured in isolation:

    X1 - T3   vs   W8 - amip_pi_base          RSNOWLIN2 alone gave trop LW +2.581
    X2 - X1   vs   U1 - T3                    DMS S=166 alone gave SO SW -2.25, trop -2.14
    X3 - X2   vs   S4 - P5                    INPPMIN 50000 gave SO SW -1.98, trop -0.10

A rung that reproduces its isolated value is additive on this base and can be reasoned
about separately.  One that does not means the base changes what the lever does, and the
configuration has to be adopted or rejected whole.

WHY THIS PAIRING WAS WORTH RETRYING.  DMS at S=166 was rejected on 2026-08-10 because it
cost the tropics 2.14 W/m2, four times the tolerance, while the tropics sit only 0.67 from
CERES period-clean.  W8 then produced +2.581 of tropical LW CRE at no resolved net energy
cost, which is the headroom the DMS lever previously had to find for itself.  Whether the
two actually compose that way is exactly what X2 measures and what nothing so far can
predict.

ONE YEAR LICENSES THE RADIATIVE METRICS ONLY.  1-yr pair thresholds, measured from the
control's interannual scatter over 1872-1915: tropical LW +-0.279, tropical SW +-0.990,
tropical net TOA +-1.656, SO SW +-1.974, SO LW +-0.584, global net TOA +-1.333.  Below
threshold is a BOUND, not a result.  No temperature claim is licensed at this length --
the seasonal T2m thresholds are 44-year numbers.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
ACC, YEAR = 3600.0, 1870
THR = {'trop LW CRE': 0.279, 'trop SW CRE': 0.990, 'trop net TOA': 1.656,
       'SO SW CRE': 1.974, 'SO LW CRE': 0.584, 'global net TOA': 1.333}
KEYS = list(THR)

LADDER = [('T3  adopted stack, S=0', 'amip_T3_dmsrev3', None),
          ('X1  + RSNOWLIN2 0.04', 'amip_X1_stack_rsnow', 'amip_T3_dmsrev3'),
          ('X2  + DMS S=166', 'amip_X2_stack_rsnow_dms166', 'amip_X1_stack_rsnow'),
          ('X3  + INPPMIN 50000', 'amip_X3_stack_rsnow_dms166_inp50k',
           'amip_X2_stack_rsnow_dms166')]
# the same lever measured on its own, for the additivity check
ISOLATED = [('RSNOWLIN2 0.04', 'amip_W8_rsnowlin204', 'amip_pi_base', 'X1  + RSNOWLIN2 0.04'),
            ('DMS S=166', 'amip_U1_dmsccn166', 'amip_T3_dmsrev3', 'X2  + DMS S=166')]

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


def row(label, d, mark=True):
    f = ''.join('*' if mark and abs(d[k]) > THR[k] else ' ' for k in KEYS)
    print(f'  {label:26s} ' + ' '.join(f'{d[k]:14.3f}' for k in KEYS) + f'  {f}')


M = {}
for lab, r, _ in LADDER:
    M[r] = metrics(r)
for _, r, c, _ in ISOLATED:
    for x in (r, c):
        if x not in M:
            M[x] = metrics(x)

print(f'year {YEAR}\n')
print(f'  {"":26s} ' + ' '.join(f'{k:>14s}' for k in KEYS))
print(f'  {"1-yr pair threshold":26s} ' + ' '.join(f'{THR[k]:14.3f}' for k in KEYS))
print('-' * 112)

print('  ABSOLUTE')
for lab, r, _ in LADDER:
    if M.get(r):
        row(lab, M[r], mark=False)

print('\n  EACH RUNG AGAINST ITS OWN CONTROL  (one variable per step)')
steps = {}
for lab, r, c in LADDER:
    if c is None or not M.get(r) or not M.get(c):
        continue
    d = {k: M[r][k] - M[c][k] for k in KEYS}
    steps[lab] = d
    row(lab, d)

print('\n  CUMULATIVE, X-rung minus T3')
base = M.get('amip_T3_dmsrev3')
for lab, r, c in LADDER[1:]:
    if base and M.get(r):
        row(lab, {k: M[r][k] - base[k] for k in KEYS})

print("\n  '*' = above its own 1-year threshold.")

# ------------------------------------------------------------- additivity
print('\nDOES EACH LEVER STILL DO WHAT IT DID ALONE?')
print('-' * 112)
for name, r, c, rung in ISOLATED:
    if not (M.get(r) and M.get(c) and rung in steps):
        print(f'  {name}: not yet available')
        continue
    iso = {k: M[r][k] - M[c][k] for k in KEYS}
    onbase = steps[rung]
    diff = {k: onbase[k] - iso[k] for k in KEYS}
    print(f'\n  {name}')
    row('   alone', iso, mark=False)
    row('   on the X base', onbase, mark=False)
    row('   difference', diff)
    big = [k for k in KEYS if abs(diff[k]) > THR[k]]
    print(f'   -> {"NOT additive on " + ", ".join(big) if big else "additive within thresholds"}')

# ------------------------------------------------------------- verdict
print('\nTHE FULL CONFIGURATION AGAINST THE CAMPAIGN TARGETS')
print('-' * 112)
for lab, r, _ in LADDER[1:]:
    if not (base and M.get(r)):
        continue
    d = {k: M[r][k] - base[k] for k in KEYS}
    tgt_lw = 100 * max(d['trop LW CRE'], 0) / 2.16
    tgt_so = 100 * max(-d['SO SW CRE'], 0) / 7.85
    trop_ok = 'within' if abs(d['trop net TOA']) <= THR['trop net TOA'] else 'OVER'
    print(f'  {lab:26s} trop LW {tgt_lw:5.0f} % of deficit | SO SW {tgt_so:5.0f} % of gap'
          f' | tropical net {d["trop net TOA"]:+6.2f} ({trop_ok})'
          f' | global {d["global net TOA"]:+6.2f}')
print('\n  Global net TOA in the T3 base is what these add to; the PI target is ~0.')
