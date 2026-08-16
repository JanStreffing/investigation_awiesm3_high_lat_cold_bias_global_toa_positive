"""Did Raupach do anything the ATMOSPHERE can see?  11J - 11I, paired, Siberian box.

WHY LOOK HERE.  Siberian soil temperature is the far end of a long causal chain, and it
carries a 1.2 K interannual noise floor: resolving the measured +0.312 K JJA soil point
estimate at 95 % would need 57 coupled years, so 11J cannot settle it even complete.
That is a statement about the DIAGNOSTIC, not about the lever.  A roughness change acts
first on momentum, and the momentum response is direct, local and much less noisy:

    z0 up  ->  |tau| up, 10 m wind down  ->  exchange coefficient up
           ->  sensible heat flux changes  ->  skin-to-air coupling changes
           ->  (eventually, weakly) soil temperature

So test the chain at its TOP, where the signal is, rather than only at the bottom where
it has been diluted into noise.  If the top of the chain is also null, Raupach genuinely
did nothing and the surface result is not merely underpowered.  If the top is strongly
resolved and the bottom is not, the lever works mechanically but its climate effect is
below what this design can detect -- a completely different conclusion, and the one that
determines whether it is worth re-testing with a better diagnostic.

NO 3-D FIELDS ARE ARCHIVED for these arms (no model-level or pressure-level output), so
"air column" here means the surface-flux and near-surface state, which is what exists.

PRE-REGISTERED SIGNS, from the offline calculation (DJF grid-box z0 0.058 -> 0.140 m,
JJA 0.190 -> 0.120 m):
    DJF:  |tau| UP, 10 m wind DOWN, |sshf| UP    (rougher)
    JJA:  |tau| DOWN, 10 m wind UP, |sshf| DOWN  (smoother -- Raupach is LOWER in summer)
The JJA sign flip is a genuine prediction and a good falsification test: noise would not
know to reverse between seasons.

ACCUMULATION TRAP.  ewss, nsss, sshf, slhf are accumulated over the output step exactly
like the radiation terms and must be divided by 3600.  Signs follow IFS convention
(downward positive for heat fluxes), so sshf is typically negative over land by day.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CTL = ('11I', f'{R092}/Tuning_test_11I_v2soil')
TRT = ('11J', f'{R092}/Tuning_test_11J_v2soil_raupach')
YEARS = list(range(1350, 1380))          # every year 11J has; 11G is not in this pair
SIB = (55.0, 75.0, 60.0, 180.0)
DJF, JJA = [12, 1, 2], [6, 7, 8]
ACC = 3600.0
# Everything IFS accumulates over the output step.  Missing one of these is a silent
# factor-3600 error, so the list is explicit rather than inferred.
FLUX = ('tsr', 'ttr', 'tsrc', 'ttrc', 'ssr', 'str', 'ssrd', 'strd',
        'ewss', 'nsss', 'sshf', 'slhf', 'e', 'tp', 'ro')


def load_year(root, var, y):
    D = f'{root}/outdata/oifs'
    div = ACC if var in FLUX else 1.0
    f = f'{D}/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(f):
        return None, None, None
    with xr.open_dataset(f, decode_times=False) as d:
        a = d[var].values / div
        lat, lon = d['lat'].values, d['lon'].values
    return (a, lat, lon) if a.shape[0] == 12 else (None, None, None)


def boxmean(f2d, lat, lon, lsm):
    la0, la1, lo0, lo1 = SIB
    ys = (lat >= la0) & (lat <= la1)
    xs = ((lon % 360) >= lo0) & ((lon % 360) <= lo1)
    sub = f2d[np.ix_(ys, xs)]
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape).copy()
    w = np.where(lsm[np.ix_(ys, xs)] >= 0.5, w, 0.0)
    return float(np.average(sub, weights=w)) if w.sum() else np.nan


def chain_for_year(root, y, lsm, lat, lon):
    g = {v: load_year(root, v, y)[0] for v in
         ('ewss', 'nsss', '10u', '10v', 'sshf', 'slhf', 'skt', '2t', 'cvh', 'lai', 'sd')}
    if g['ewss'] is None or g['10u'] is None:
        return None
    out = {}
    for tag, mons in (('DJF', DJF), ('JJA', JJA)):
        m = [x - 1 for x in mons]
        tau = np.sqrt(g['ewss'][m] ** 2 + g['nsss'][m] ** 2).mean(axis=0)
        wind = np.sqrt(g['10u'][m] ** 2 + g['10v'][m] ** 2).mean(axis=0)
        out[f'{tag} |tau| [N/m2]'] = boxmean(tau, lat, lon, lsm)
        out[f'{tag} 10m wind [m/s]'] = boxmean(wind, lat, lon, lsm)
        if g['sshf'] is not None:
            out[f'{tag} sshf [W/m2]'] = boxmean(g['sshf'][m].mean(axis=0), lat, lon, lsm)
        if g['slhf'] is not None:
            out[f'{tag} slhf [W/m2]'] = boxmean(g['slhf'][m].mean(axis=0), lat, lon, lsm)
        if g['skt'] is not None and g['2t'] is not None:
            out[f'{tag} skt-2t [K]'] = boxmean((g['skt'] - g['2t'])[m].mean(axis=0),
                                               lat, lon, lsm)
        if g['cvh'] is not None:
            out[f'{tag} cvh [-]'] = boxmean(g['cvh'][m].mean(axis=0), lat, lon, lsm)
        if g['sd'] is not None:
            out[f'{tag} snow depth [m]'] = boxmean(g['sd'][m].mean(axis=0), lat, lon, lsm)
    return out


def series(root):
    lsm, lat, lon = load_year(root, 'lsm', YEARS[0])
    if lsm is not None and lsm.ndim == 3:
        lsm = lsm[0]
    out = {}
    for y in YEARS:
        r = chain_for_year(root, y, lsm, lat, lon)
        if r:
            out[y] = r
    return out


print(__doc__)
print('=' * 104)
sc, st = series(CTL[1]), series(TRT[1])
yrs = [y for y in YEARS if y in sc and y in st]
print(f'\npaired over {len(yrs)} years, {yrs[0]}-{yrs[-1]}, Siberian land '
      f'{SIB[0]:.0f}-{SIB[1]:.0f}N {SIB[2]:.0f}-{SIB[3]:.0f}E\n')

try:
    from scipy import stats
    HAVE = True
except ImportError:
    HAVE = False

keys = [k for k in sc[yrs[0]] if k in st[yrs[0]]]
order = [k for tag in ('DJF', 'JJA') for k in keys if k.startswith(tag)]
print(f'  {"metric":24s} {CTL[0]:>10s} {TRT[0]:>10s} {"diff":>9s} {"%":>7s} '
      f'{"t":>7s} {"p":>8s}   {"95% CI":>19s}')
res = {}
for k in order:
    a = np.array([sc[y][k] for y in yrs])
    b = np.array([st[y][k] for y in yrs])
    d = b - a
    n = d.size
    se = d.std(ddof=1) / np.sqrt(n)
    t = d.mean() / se if se > 0 else np.nan
    if HAVE:
        p = 2 * (1 - stats.t.cdf(abs(t), n - 1))
        tc = stats.t.ppf(0.975, n - 1)
    else:
        p, tc = np.nan, 2.045
    pct = 100.0 * d.mean() / abs(a.mean()) if a.mean() != 0 else np.nan
    res[k] = (d.mean(), t, p)
    star = '*' if (p == p and p < 0.05) else ' '
    if k.startswith('JJA') and order.index(k) and not order[order.index(k) - 1].startswith('JJA'):
        print()
    print(f'  {k:24s} {a.mean():10.4f} {b.mean():10.4f} {d.mean():+9.4f}{star}'
          f' {pct:+6.1f}% {t:7.2f} {p:8.4f}   [{d.mean()-tc*se:+8.4f},{d.mean()+tc*se:+8.4f}]')

if HAVE:
    m = len(order)
    print(f'\n  Bonferroni over {m} tests: alpha = {0.05/m:.4f}')
    for k in order:
        if res[k][2] < 0.05 / m:
            print(f'    SURVIVES  {k:24s} d={res[k][0]:+.4f}  p={res[k][2]:.5f}')

print('\n  PRE-REGISTERED SIGN TEST (noise cannot know to flip between seasons)')
for a_key, b_key, want in (('DJF |tau| [N/m2]', 'JJA |tau| [N/m2]', 'DJF up, JJA down'),
                           ('DJF 10m wind [m/s]', 'JJA 10m wind [m/s]',
                            'DJF down, JJA up')):
    if a_key in res and b_key in res:
        da, db = res[a_key][0], res[b_key][0]
        ok = (da > 0 > db) if 'tau' in a_key else (da < 0 < db)
        print(f'    {a_key.split()[1]:12s} DJF {da:+.4f}  JJA {db:+.4f}   '
              f'expected {want}: {"MATCHES" if ok else "does not match"}')
print()
