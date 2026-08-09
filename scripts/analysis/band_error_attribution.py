"""Where is each band's radiative error, and which lever moves which band?

WHY THIS EXISTS.  Round 23 was designed to close the Southern Ocean cloud-AMOUNT
deficit, and every candidate in it (less cloud erosion, delayed maritime warm rain,
a higher INP gate) adds low cloud and therefore REFLECTS more.  That cools.  But the
tropics are already too cold and already absorb too little, so a lever that is right
for the Southern Ocean may be unadoptable because of what it does 40 degrees away.
The campaign has hit this before -- A1a was rejected exactly this way.

So before running anything else, two questions need answering properly:

  1. WHAT KIND of error does each band have?  A band that absorbs too much because its
     CLEAR-SKY albedo is too low needs a different fix from one whose CLOUD does not
     reflect enough, and both look identical in net TOA.  The four-way split
     (SW clear, SW cloud, LW clear, LW cloud) is the minimum needed to tell them apart.

  2. WHICH LEVERS WARM THE TROPICS?  If an SO lever costs the tropics, it can only be
     adopted alongside something that pays that cost back.  The campaign has 50 runs and
     has never asked which of them move the tropics in the useful direction.  This is
     the compensating-bias question: two errors of opposite sign can both be real, and
     fixing one alone makes the total worse.

METHOD.  amip_presentday 1990-2014 against CERES EBAF 07/2005-06/2015, period-clean, so
the epoch offset that makes the PI arm misleading (it turns the tropical error from
-0.67 into -2.51) does not apply.  Band means are exact on each source's own grid, so
nothing is remapped.  The lever table necessarily uses the PI arms, where the levers
were run, but lever DELTAS are epoch-insensitive because both sides share the forcing.

SIGN CONVENTION throughout: positive = the model gains energy relative to CERES, i.e.
warms.  So a positive SW-cloud entry means the model's cloud reflects too LITTLE.
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
PD, PDY = 'amip_presentday', (1990, 2014)
LEVY = (1872, 1915)

BANDS = [('90S-65S', -90, -65), ('SO 65S-45S', -65, -45), ('45S-30S', -45, -30),
         ('tropics', -30, 30), ('30N-45N', 30, 45), ('45N-65N', 45, 65),
         ('65N-90N', 65, 90)]

LEVERS = [('A1a ovlliqice0.10', 'amip_A1_overlap01'),
          ('A1b ovlliqice0.35', 'amip_A1_overlap035'),
          ('A1c depliqdep1500', 'amip_A1c_depliqdepth1500'),
          ('B2 clddiffconvi25', 'amip_B2_clddiffconvi25'),
          ('B3 clddiff 1.5e-5', 'amip_B3_clddiff15e6'),
          ('B4 entshalp3', 'amip_B4_entshalp3'),
          ('B5 capdcycl0', 'amip_B5_capdcycl0'),
          ('B6 lcritsnow1e-5', 'amip_B6_lcritsnow1e5'),
          ('B7 rvice 0.22', 'amip_B7_rvice022'),
          ('D1 capdcycl4', 'amip_D1_capdcycl4'),
          ('D2a inpsea 0.2', 'amip_D2a_inpsea02'),
          ('D2b inp+p700', 'amip_D2b_inpsea02_p700'),
          ('F4 rsmin1000', 'amip_F4_rsmin1000'),
          ('G4 tundra225', 'amip_G4_tundra'),
          ('K1 landalb', 'amip_K1_landalb')]

print(__doc__)
print('=' * 108)


def load(run, y0, y1, vars_):
    D = f'{RT}/{run}/outdata/oifs'
    out, lat = {}, None
    for v in vars_:
        acc = []
        for y in range(y0, y1 + 1):
            f = f'{D}/atm_remapped_1m_{v}_1m_{y}-{y}.nc'
            if not os.path.exists(f):
                continue
            with xr.open_dataset(f, decode_times=False) as d:
                a = d[v].values
                if lat is None:
                    lat = d['lat'].values
            acc.append(a.mean(axis=0) / (1.0 if v in ('tcc', 'lcc', 'mcc', 'hcc') else ACC))
        if not acc:
            return None, None
        out[v] = np.mean(acc, axis=0)
    return out, lat


def bm(f2d, lat, a, b):
    s = (lat >= a) & (lat < b)
    return float(np.average(f2d[s, :].mean(axis=1), weights=np.cos(np.deg2rad(lat[s]))))


def afrac(lat, a, b):
    w = np.cos(np.deg2rad(lat))
    return float(w[(lat >= a) & (lat < b)].sum() / w.sum())


cds = xr.open_dataset(CERESF)
clat = cds['lat'].values


def cb(v, a, b):
    s = (clat >= a) & (clat < b)
    x = cds[v].values.mean(axis=0)
    return float(np.average(x[s, :].mean(axis=1), weights=np.cos(np.deg2rad(clat[s]))))


V = ['tsr', 'ttr', 'tsrc', 'ttrc', 'tcc', 'hcc', 'lcc']
m, lat = load(PD, *PDY, V)

print('1. FOUR-WAY ATTRIBUTION OF EACH BAND, period-clean.  '
      'Positive = model gains energy (warms) vs CERES.')
print('-' * 108)
print(f'  {"band":12s} {"net":>7s} | {"SWclr":>7s} {"SWcld":>7s} {"LWclr":>7s} {"LWcld":>7s} | '
      f'{"dCLD":>6s} {"dHIGH":>6s} | {"contrib":>8s}   dominant term')
tot = 0.0
rows = []
for nm, a, b in BANDS:
    # model
    mnet = bm(m['tsr'] + m['ttr'], lat, a, b)
    mswc = bm(m['tsr'] - m['tsrc'], lat, a, b)       # SW CRE
    mlwc = bm(m['ttr'] - m['ttrc'], lat, a, b)       # LW CRE
    mswclr = bm(m['tsrc'], lat, a, b)
    mlwclr = bm(m['ttrc'], lat, a, b)
    mcld = bm(m['tcc'], lat, a, b) * 100
    mhigh = bm(m['hcc'], lat, a, b) * 100
    # CERES: clear-sky net = net - CRE(sw) - CRE(lw); use total-region clear-sky
    onet = cb('toa_net_all_clim', a, b)
    oswc = cb('toa_cre_sw_clim', a, b)
    olwc = cb('toa_cre_lw_clim', a, b)
    olwup_clr = cb('toa_lw_clr_t_clim', a, b)
    oswup_clr = cb('toa_sw_clr_t_clim', a, b)
    ocld = cb('cldarea_total_daynight_clim', a, b)
    # clear-sky absorbed SW = net_clr + OLR_clr ; net_clr = net - swcre - lwcre
    onet_clr = onet - oswc - olwc
    oswclr = onet_clr + olwup_clr
    olwclr = -olwup_clr

    d_swclr = mswclr - oswclr
    d_swcld = mswc - oswc
    d_lwclr = mlwclr - olwclr
    d_lwcld = mlwc - olwc
    dnet = mnet - onet
    fr = afrac(lat, a, b)
    tot += dnet * fr
    terms = {'SWclr': d_swclr, 'SWcld': d_swcld, 'LWclr': d_lwclr, 'LWcld': d_lwcld}
    dom = max(terms, key=lambda k: abs(terms[k]))
    rows.append((nm, dnet, terms, dnet * fr, mcld - ocld, mhigh))
    print(f'  {nm:12s} {dnet:+7.2f} | {d_swclr:+7.2f} {d_swcld:+7.2f} {d_lwclr:+7.2f} '
          f'{d_lwcld:+7.2f} | {mcld-ocld:+6.2f} {mhigh:6.1f} | {dnet*fr:+8.3f}   {dom}')
print(f'  {"GLOBAL":12s} {"":>7s} |{"":>32s} | {"":>13s} | {tot:+8.3f}')

print('\n  Read the two SW columns against each other: SWclr is surface/aerosol albedo,')
print('  SWcld is how much the cloud reflects.  LWcld is high-cloud greenhouse trapping.')

trop = [r for r in rows if r[0] == 'tropics'][0]
print(f"\n2. THE TROPICS SPECIFICALLY: net {trop[1]:+.2f} W/m2")
print('-' * 108)
for k, v in trop[2].items():
    print(f'    {k:8s} {v:+7.2f}')
worst = min(trop[2], key=lambda k: trop[2][k])
print(f'  Largest COOLING term: {worst} at {trop[2][worst]:+.2f} W/m2.')
if worst == 'LWcld':
    print('  => The tropics lose energy because their high cloud traps too little longwave.')
    print('     With SST prescribed and correct, that is an anvil/cirrus deficit -- too')
    print('     little ice cloud, too low, or too thin.  A lever that puts ice cloud BACK')
    print('     warms the tropics WITHOUT touching the Southern Ocean low-cloud problem,')
    print('     because it acts in a different phase, a different altitude and a different')
    print('     regime.  That is the disjointness that made F4+D2b superpose.')

print('\n3. WHICH EXISTING LEVERS WARM THE TROPICS?  '
      'Delta vs control, PI arms, 44 yr.')
print('-' * 108)
base, blat = load('amip_pi_base', *LEVY, ['tsr', 'ttr', 'tsrc', 'ttrc', 'hcc'])
bnet_t = bm(base['tsr'] + base['ttr'], blat, -30, 30)
blw_t = bm(base['ttr'] - base['ttrc'], blat, -30, 30)
bsw_t = bm(base['tsr'] - base['tsrc'], blat, -30, 30)
bhi_t = bm(base['hcc'], blat, -30, 30) * 100
bnet_so = bm(base['tsr'] + base['ttr'], blat, -65, -45)
bsw_so = bm(base['tsr'] - base['tsrc'], blat, -65, -45)

print(f'  {"lever":19s} {"dTROP":>7s} {"dLWcre":>7s} {"dSWcre":>7s} {"dHIGH":>6s} | '
      f'{"dSO":>7s} {"dSOcre":>7s} | verdict')
cands = []
for nm, run in LEVERS:
    d, dl = load(run, *LEVY, ['tsr', 'ttr', 'tsrc', 'ttrc', 'hcc'])
    if d is None:
        continue
    dt = bm(d['tsr'] + d['ttr'], dl, -30, 30) - bnet_t
    dlw = bm(d['ttr'] - d['ttrc'], dl, -30, 30) - blw_t
    dsw = bm(d['tsr'] - d['tsrc'], dl, -30, 30) - bsw_t
    dhi = bm(d['hcc'], dl, -30, 30) * 100 - bhi_t
    dso = bm(d['tsr'] + d['ttr'], dl, -65, -45) - bnet_so
    dsocre = bm(d['tsr'] - d['tsrc'], dl, -65, -45) - bsw_so
    v = ''
    if dt > 0.2 and dso < -0.2:
        v = '*** WARMS TROPICS *AND* DARKENS SO'
    elif dt > 0.2:
        v = 'warms tropics'
    elif dso < -1.0:
        v = 'good for SO, costs tropics'
    cands.append((nm, dt, dlw, dsw, dhi, dso, dsocre, v))
    print(f'  {nm:19s} {dt:+7.2f} {dlw:+7.2f} {dsw:+7.2f} {dhi:+6.2f} | '
          f'{dso:+7.2f} {dsocre:+7.2f} | {v}')

print('\n  dSO negative = the Southern Ocean absorbs less, which is the direction needed.')
print('  dTROP positive = the tropics absorb more, which is also the direction needed.')
warm = [c for c in cands if c[1] > 0.2]
if warm:
    best = max(warm, key=lambda c: c[1])
    print(f'\n  Tropics-warming levers found: {", ".join(c[0] for c in warm)}')
    print(f'  Largest: {best[0]} at {best[1]:+.2f} W/m2, via dLWcre {best[2]:+.2f} / '
          f'dSWcre {best[3]:+.2f}, high cloud {best[4]:+.2f} pp.')
else:
    print('\n  *** NO LEVER IN 50 RUNS WARMS THE TROPICS by more than 0.2 W/m2.')
    print('      Every tropical response on record is neutral or cooling, so the pairing')
    print('      the Southern Ocean needs does not exist yet and has to be built.')
cds.close()
