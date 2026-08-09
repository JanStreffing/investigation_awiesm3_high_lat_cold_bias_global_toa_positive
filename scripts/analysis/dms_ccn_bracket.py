"""The marine-biogenic-CCN bracket: how much Southern Ocean can be bought, and at what
cost to the tropics?

THE LEVER.  radpar.F90 now adds a DMS-derived term to PCCNO over ocean:

    Sc  Saltzman et al. (1993);  k  Nightingale et al. (2000);  F = k * Cw
    PCCNO += ECE_DMS_CCN_SENS * F_S        [F_S in mg S m-2 d-1]

PCCNO flows through the Martin et al. (1994) quadratic to droplet number and effective
radius, i.e. the Twomey effect.  The input is DMS-Rev3 (Hulswar et al. 2022), ingested
through ocp-tool and verified in T3.

WHY A BRACKET AND NOT A VALUE.  Woodhouse et al. (2010, ACP 10, 7545) measure
dCCN/dFlux_DMS anywhere from -43 to +166 cm-3 per (mg S m-2 d-1) depending on season and
hemisphere.  That range is far too wide to pick a number from honestly, so the
coefficient is a namelist knob and the range is SAMPLED: T3 (S=0), U2 (43), U3 (90),
U1 (166).  Four points, one control, all sharing identical DMS input and identical winds.

WHAT THE CURVE IS FOR.  U1 already showed the lever is real (SO SW CRE -2.25 W/m2, above
the 1-year detection threshold of 1.97) but that it spends the tropics at four times the
pre-registered guardrail.  The question is therefore not "does it work" but "is there any
S in the published range that buys useful Southern Ocean without wrecking the tropics".
With four points that is READ OFF a curve rather than guessed at with a fifth run.

The Twomey response goes as d ln(Nd), so it must be SUB-LINEAR in S: doubling the
coefficient buys less than double.  That is checked here, not assumed -- and it matters
in the unfavourable direction, because it means halving S from a failing value removes
MORE than half the SO benefit while removing only half the tropical cost... or the
reverse.  Which one it is, is exactly what a fitted curve settles.

SEASONALITY IS THE WHOLE ARGUMENT.  U1's annual means look like a wash (SO -2.25 against
tropics -2.14) but that hides the structure: the SO response is concentrated in DJF
(-7.25) because that is when the bloom is, while the tropical cost is flat year-round.
The SO deficit is also worst in DJF (+10.86 W/m2).  So the seasonal table below is the
decisive output, not the annual one -- and the selectivity comes from the DMS field
itself, with no hemispheric parameter invented to produce it.

THRESHOLDS, fixed in advance from interannual scatter (round 22/23):
    SO SW CRE, 1-year pair   +-1.97 W/m2      <- can this run see the effect at all
    SO SW CRE, 44-year pair  +-0.30 W/m2
    tropical guardrail       +-0.50 W/m2      <- the tropics are only -0.67 from CERES
                                                 period-clean; more than this and the
                                                 lever is a global knob in disguise

CAVEAT, stated before the numbers.  One year resolves a cloud lever only above ~2 W/m2
(the one-year screening rule, round 24).  U1 clears it; U2 at the low end may not, in
which case its SO number BOUNDS the effect rather than resolving it.  Flagged per row.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

ACC = 3600.0                      # IFS radiative fluxes are accumulated J/m2 per hour
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
YEAR = 1870                       # the single screening year all four runs share
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'

BANDS = [('SO 65-45S', -65.0, -45.0),
         ('tropics 30S-30N', -30.0, 30.0),
         ('global', -90.0, 90.0)]
SEAS = {'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8], 'SON': [9, 10, 11]}

# S = 0 is the control: identical DMS-Rev3 input, identical winds, lever off.
ARMS = [(0.0, 'amip_T3_dmsrev3'), (43.0, 'amip_U2_dmsccn43'),
        (90.0, 'amip_U3_dmsccn90'), (166.0, 'amip_U1_dmsccn166')]

THR_SO_1Y, THR_SO_44Y, GUARD_TROP = 1.97, 0.30, 0.50
SO_DJF_DEFICIT = 10.86            # model minus CERES, DJF SO SW CRE, period-clean

print(__doc__)
print('=' * 100)


def swcre(run):
    """(12, nlat, nlon) monthly SW cloud radiative effect [W/m2], and lat."""
    D = f'{RT}/{run}/outdata/oifs'
    out = {}
    for var in ('tsr', 'tsrc'):
        f = f'{D}/atm_remapped_1m_{var}_1m_{YEAR}-{YEAR}.nc'
        with xr.open_dataset(f, decode_times=False) as d:
            out[var] = d[var].values / ACC
            lat = d['lat'].values
    return out['tsr'] - out['tsrc'], lat


def band(a3d, lat, lo, hi, months=None):
    """Cosine-weighted band mean over the given months (all 12 if None)."""
    a = a3d if months is None else a3d[[m - 1 for m in months]]
    sel = (lat >= lo) & (lat < hi)
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(a.mean(axis=0)[sel, :].mean(axis=1), weights=w))


# ------------------------------------------------------------------ load the four arms
arms = []
for s, run in ARMS:
    try:
        a, lat = swcre(run)
    except FileNotFoundError:
        print(f'  MISSING: {run} -- {YEAR} not on disk, skipped')
        continue
    arms.append((s, run, a, lat))
if len(arms) < 2:
    raise SystemExit('need the control plus at least one arm')
ctrl_s, ctrl_run, ctrl, lat = arms[0]
assert ctrl_s == 0.0, 'first arm must be the S=0 control'

print(f'control {ctrl_run} (S=0), year {YEAR}; arms: '
      f'{", ".join(f"S={s:.0f}" for s, _, _, _ in arms[1:])}\n')

# ------------------------------------------------------------------ 1. annual response
print('1. ANNUAL SW CRE RESPONSE, each arm minus the S=0 control [W/m2]')
print('-' * 100)
hdr = f'  {"S":>5s} ' + ' '.join(f'{n:>18s}' for n, _, _ in BANDS) + '   verdict'
print(hdr)
print(f'  {"0":>5s} ' + ' '.join(f'{band(ctrl, lat, a, b):18.2f}' for _, a, b in BANDS)
      + '   (control, absolute)')
annual = {}
for s, run, a, _ in arms[1:]:
    d = {n: band(a, lat, lo, hi) - band(ctrl, lat, lo, hi) for n, lo, hi in BANDS}
    annual[s] = d
    so, tr = d['SO 65-45S'], d['tropics 30S-30N']
    seen = 'resolved' if abs(so) > THR_SO_1Y else 'BELOW 1-yr floor (bound only)'
    guard = 'within guardrail' if abs(tr) <= GUARD_TROP else \
            f'GUARDRAIL x{abs(tr) / GUARD_TROP:.1f}'
    print(f'  {s:5.0f} ' + ' '.join(f'{d[n]:18.2f}' for n, _, _ in BANDS)
          + f'   SO {seen}; tropics {guard}')

# ------------------------------------------------------------------ 2. seasonality
print('\n2. SEASONAL RESPONSE -- the SO gain is concentrated, the tropical cost is not')
print('-' * 100)
print(f'  {"S":>5s} {"season":>7s} {"SO dCRE":>10s} {"trop dCRE":>10s} {"ratio":>8s}'
      f'  {"% of DJF SO deficit":>21s}')
season_tab = {}
for s, run, a, _ in arms[1:]:
    season_tab[s] = {}
    for sn, mm in SEAS.items():
        so = band(a, lat, -65, -45, mm) - band(ctrl, lat, -65, -45, mm)
        tr = band(a, lat, -30, 30, mm) - band(ctrl, lat, -30, 30, mm)
        season_tab[s][sn] = (so, tr)
        ratio = abs(so / tr) if abs(tr) > 1e-9 else float('inf')
        closed = f'{100 * abs(so) / SO_DJF_DEFICIT:19.0f} %' if sn == 'DJF' else ' ' * 21
        print(f'  {s:5.0f} {sn:>7s} {so:10.2f} {tr:10.2f} {ratio:8.2f}  {closed}')
    print()

# ------------------------------------------------------------------ 3. the response curve
print('3. THE RESPONSE CURVE -- is it linear in S, and where does the guardrail bind?')
print('-' * 100)
S = np.array([s for s, _, _, _ in arms[1:]])
if len(S) >= 2:
    for label, key, sub in (('SO annual', 'SO 65-45S', None),
                            ('SO DJF', None, 'DJF'),
                            ('tropics annual', 'tropics 30S-30N', None),
                            ('tropics DJF', None, 'DJF')):
        if sub is None:
            y = np.array([annual[s][key] for s in S])
        else:
            idx = 0 if label.startswith('SO') else 1
            y = np.array([season_tab[s][sub][idx] for s in S])
        # per-unit-S response at each sampled point: flat => linear, falling => saturating
        per = y / S
        shape = ('linear in S' if per.max() - per.min() < 0.1 * abs(per.mean())
                 else 'SUB-linear (saturating)' if abs(per[-1]) < abs(per[0])
                 else 'SUPER-linear')
        print(f'  {label:16s} ' + ' '.join(f'S={s:.0f}: {v:7.2f}' for s, v in zip(S, y))
              + f'   | per unit S: ' + ' '.join(f'{p:6.3f}' for p in per)
              + f'   -> {shape}')

    # Is the SO response even ordered in S?  If it is not, these are noise, not a curve,
    # and no coefficient can be read off them.  This check has to come first: np.interp
    # would happily CLAMP to the sampled range and report a confident wrong answer.
    yso = np.array([annual[s]['SO 65-45S'] for s in S])
    yso_djf = np.array([season_tab[s]['DJF'][0] for s in S])
    ytr = np.array([abs(annual[s]['tropics 30S-30N']) for s in S])
    mono_so = bool(np.all(np.diff(yso) <= 0) or np.all(np.diff(yso) >= 0))
    mono_tr = bool(np.all(np.diff(ytr) >= 0) or np.all(np.diff(ytr) <= 0))
    resolved = [s for s in S if abs(annual[s]['SO 65-45S']) > THR_SO_1Y]
    print(f'\n  SO response monotonic in S?      {"yes" if mono_so else "NO"}'
          f'      arms resolved at 1 yr: {[int(s) for s in resolved] or "none but S=166"}')
    print(f'  tropical response monotonic in S? {"yes" if mono_tr else "NO"}')
    if not mono_so:
        print('  => the SO points are NOT a curve.  Every arm below the 1-year floor is\n'
              '     internal variability, and the non-monotonicity is the proof of it.\n'
              '     No coefficient can be interpolated from them.')

    # The tropical response is the one that can be extrapolated, if it is ordered.  Do it
    # through the origin (S=0 must give zero response by construction), and say plainly
    # when the answer falls outside the sampled range instead of clamping into it.
    if mono_tr:
        slope = float(np.polyfit(np.r_[0.0, S], np.r_[0.0, ytr], 1)[0])
        s_guard = GUARD_TROP / slope
        where = ('inside' if 43 <= s_guard <= 166 else
                 'BELOW the whole published range' if s_guard < 43 else 'above it')
        print(f'\n  tropical cost is {slope:.4f} W/m2 per unit S (fit through the origin)')
        print(f'  => guardrail |tropics| <= {GUARD_TROP} binds at S ~ {s_guard:.0f}, '
              f'which is {where} (Woodhouse 43..166)')
        if s_guard < 43:
            print(f'  At the lowest published value, S=43, the tropics already move '
                  f'{annual[43.0]["tropics 30S-30N"]:+.2f} W/m2\n'
                  f'  while the SO DJF gain, {yso_djf[0]:+.2f}, is below what one year '
                  f'can resolve.')

print('\n' + '=' * 100)
print("""VERDICT LOGIC (fixed before the runs):
  tropics within +-0.5 AND |SO DJF| > 2  -> adoptable; promote to 44 yr
  tropics within +-0.5 AND |SO DJF| < 2  -> real but too small to matter; retire
  tropics beyond +-0.5 at every S >= 43  -> no value in the published range works""")
