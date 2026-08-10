"""Which lever actually cools the Southern Ocean?  Ranked on NET, not shortwave.

WHY THIS EXISTS.  This campaign -- and I did it repeatedly on 2026-08-10 -- scores
Southern Ocean levers on SW CRE alone.  That is the wrong scorecard.  The coupled model's
complaint is a WARM SST bias, which responds to the NET energy the cloud lets through:

    d(SO net cloud forcing)  =  d(SW CRE)  +  d(LW CRE)

A lever that reflects more sunlight AND traps more outgoing longwave has bought much less
than its shortwave number suggests, and one that does both equally has bought nothing.
RSNOWLIN2 is the case in point: +2.6 to +3.6 W/m2 of SO LW CRE, which cancels most of its
shortwave gain and in some configurations flips the net POSITIVE, i.e. it warms the
Southern Ocean while looking like an improvement on the usual metric.

WHAT IS RANKED.  Every 44-year arm in the archive, from the SO SW and SO LW columns
already computed by tropical_lw_scan.py, plus the one-year screens from rounds 26-27
scored against their own controls.  The two groups are NEVER mixed in one ranking: their
thresholds differ by a factor of six.

    44-year pair:  SO SW +-0.298   SO LW +-0.088   SO net TOA +-0.296
    1-year pair:   SO SW +-1.974   SO LW +-0.584   SO net TOA +-1.967

CAVEAT ON THE PROXY.  SW CRE + LW CRE is the CLOUD contribution to the net, not the full
net TOA change -- it omits clear-sky changes.  For cloud microphysics levers that is
negligible; for anything touching aerosol or surface albedo it is not, and those rows are
flagged.  The one-year block computes true net TOA directly and is not a proxy.

READING IT.  Negative net = less energy absorbed = cooling the Southern Ocean = the
direction we want.  The band is 7.85 W/m2 short of CERES in SW CRE, but what has to be
closed in the NET is smaller, because the model's SO LW CRE is not itself badly wrong.
Costs are printed alongside: a lever is only useful if the tropics and the global budget
survive it.
"""
import os, re, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BASE = ('/work/ab0246/a270092/postprocessing/'
        'investigation_awiesm3_high_lat_cold_bias_global_toa_positive')
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
ACC, Y1YR = 3600.0, 1870
SCAN = f'{BASE}/data/tropical_lw_scan.out'
# levers whose clear-sky response is NOT negligible, so the CRE proxy understates them
AEROSOL_OR_SURFACE = ('K1', 'K2', 'G2', 'G3', 'G4', 'F1', 'F2', 'F3', 'F4', 'F5',
                      'H1', 'H2', 'I1', 'I2', 'I3', 'J1', 'J2', 'L1', 'L2', 'N2',
                      'O1', 'O2', 'P3', 'P4', 'P5', 'P6', 'expA', 'piCTRL')
THR44 = {'SO SW': 0.298, 'SO LW': 0.088, 'SO net': 0.296}
THR1 = {'SO SW': 1.974, 'SO LW': 0.584, 'SO net': 1.967, 'trop net': 1.656,
        'glob net': 1.333}

print(__doc__)
print('=' * 104)

# ------------------------------------------------------------ A. the 44-year archive
print('A. THE 44-YEAR ARCHIVE, ranked by NET cloud forcing over 45-65S')
print('-' * 104)
rows = []
if os.path.exists(SCAN):
    started = False
    for line in open(SCAN):
        if line.startswith('  run '):
            started = True
            continue
        if not started:
            continue
        if line.startswith('\n') or line.strip().startswith('not scored'):
            break
        m = re.match(r'\s{2}(\S.*?)\s{2,}(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
                     r'\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', line)
        if m:
            name = m.group(1).strip()
            tlw, tsw, sosw, solw, gnet = (float(m.group(i)) for i in range(2, 7))
            rows.append((name, sosw, solw, sosw + solw, gnet))
if not rows:
    print('  tropical_lw_scan.out not available -- run tropical_lw_scan.py first')
else:
    rows.sort(key=lambda r: r[3])
    print(f'  {"run":22s} {"SO SW":>9s} {"SO LW":>9s} {"SO NET":>9s} {"global net":>11s}'
          f'   flags')
    for name, sw, lw, net, g in rows[:14]:
        f = []
        if abs(net) > THR44['SO net']:
            f.append('net resolved')
        if lw > THR44['SO LW'] and sw < -THR44['SO SW']:
            f.append(f'LW cancels {100*min(1,lw/abs(sw)):.0f}%')
        if any(name.startswith(a) for a in AEROSOL_OR_SURFACE):
            f.append('CRE proxy understates')
        print(f'  {name:22s} {sw:9.3f} {lw:9.3f} {net:9.3f} {g:11.3f}   {"; ".join(f)}')
    print(f'\n  ({len(rows)} arms scored; showing the 14 that cool the SO most in NET)')
    worst = [r for r in rows if r[3] > THR44['SO net']]
    if worst:
        print(f'  WARMING the SO in net, above threshold: '
              f'{", ".join(f"{n} ({v:+.2f})" for n, _, _, v, _ in worst[-4:])}')

# ------------------------------------------------------------ B. the one-year screens
print('\nB. ROUNDS 26-27 ONE-YEAR SCREENS, true net TOA, each against its own control')
print('-' * 104)
PAIRS = [('W8 RSNOWLIN2 .04', 'amip_W8_rsnowlin204', 'amip_pi_base'),
         ('W4 DETRPEN 1.32e-4', 'amip_W4_detrpen132', 'amip_pi_base'),
         ('W6 RVICE 0.18', 'amip_W6_rvice018', 'amip_pi_base'),
         ('W9 DARS2 stack', 'amip_W9_dars2stack', 'amip_pi_base'),
         ('X1 stack+rsnow', 'amip_X1_stack_rsnow', 'amip_T3_dmsrev3'),
         ('X4 +inp50k', 'amip_X4_stack_rsnow_inp50k', 'amip_T3_dmsrev3'),
         ('X3 +dms+inp50k', 'amip_X3_stack_rsnow_dms166_inp50k', 'amip_T3_dmsrev3'),
         ('Y1 +ovl 0.35', 'amip_Y1_ovl035', 'amip_T3_dmsrev3'),
         ('Y2 +ovl 0.10', 'amip_Y2_ovl01', 'amip_T3_dmsrev3'),
         ('  DMS alone (U1)', 'amip_U1_dmsccn166', 'amip_T3_dmsrev3')]


def met(run):
    d = {}
    for v in ('tsr', 'tsrc', 'ttr', 'ttrc'):
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{v}_1m_{Y1YR}-{Y1YR}.nc'
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
    return {'SO SW': b(sw, -65, -45), 'SO LW': b(lw, -65, -45),
            'SO net': b(net, -65, -45), 'trop net': b(net, -30, 30),
            'glob net': b(net, -90, 90)}


cache = {}
out = []
for lab, r, c in PAIRS:
    for x in (r, c):
        if x not in cache:
            cache[x] = met(x)
    if not cache[r] or not cache[c]:
        continue
    d = {k: cache[r][k] - cache[c][k] for k in THR1}
    out.append((lab, d))
out.sort(key=lambda t: t[1]['SO net'])
print(f'  {"arm":22s} {"SO SW":>9s} {"SO LW":>9s} {"SO NET":>9s} {"trop net":>10s}'
      f' {"glob net":>10s}   resolved')
for lab, d in out:
    r = ''.join('*' if abs(d[k]) > THR1[k] else ' '
                for k in ('SO SW', 'SO LW', 'SO net', 'trop net', 'glob net'))
    print(f'  {lab:22s} {d["SO SW"]:9.3f} {d["SO LW"]:9.3f} {d["SO net"]:9.3f} '
          f'{d["trop net"]:10.3f} {d["glob net"]:10.3f}   {r}')
print('\n  columns of the resolved flag: SO SW | SO LW | SO net | trop net | glob net')
print('  1-yr thresholds: 1.974 / 0.584 / 1.967 / 1.656 / 1.333')
