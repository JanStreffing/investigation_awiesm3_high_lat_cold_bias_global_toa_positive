"""DMS coupled: 11F against 11E, the 20 years that exist.

WHY A CANCELLED RUN IS BEING SCORED.  11F was killed at 20 of 50 years because S4 had
superseded its lever within hours of submission: nearly the same annual Southern Ocean SW
CRE forcing (S4 -1.97 against DMS -2.25) at a tropical cost of 0.10 against DMS's 2.14,
namelist-only, and already run.  The generic question it was kept for -- what ~2 W/m^2 of
SO brightening does once SST and sea ice can RESPOND -- was judged not worth 30 more years.

That judgement rested on S4 closing the Southern Ocean.  It does not.  Scored coupled in
round 28, S4 removes only 27 % of 11E's CERES gap (+5.36 -> +3.89), and its coupled
response is 0.28-0.31x its AMIP response, not the 1.40x predicted from the DMS pair.  The
question is live again and 20 paired years are already on disk.

THE PAIR IS CLEAN, checked rather than assumed -- the S4 comparison lost a decade to
exactly this.  Verified 2026-08-18, both legs 1350-59 and 1360-69:

    staged slt_TCO95.nc organic (type 6) cells   11E 547        11F 547        matched
    LPJ-GUESS binary md5                         11E 8c5ab467   11F 8c5ab467   matched

So 11F is 11E with ECE_DMS_CCN_SENS 0 -> 166 and nothing else.

WHAT THIS CAN AND CANNOT SETTLE.  It CANNOT rehabilitate DMS as a tuning knob.  The
tropical cost is 0.0129 W/m^2 per unit S, so the +-0.5 tolerance binds at S~39, below
Woodhouse et al.'s published floor of 43; and scaling S cannot separate the bands because
tropical background CCN is LOWER than the SO's.  That is forcing selectivity, and free
SSTs do not touch it.

What it CAN settle is the screening rule for every future SO lever: does a lever with a
large SO SURFACE response amplify or damp under coupling?  S4 damped at 0.28x while
growing 0.595 pp of sea ice (not resolved).  If DMS amplifies here while growing ~3.2 pp,
then "amplification is a property of the lever, carried by sea-ice albedo" rests on two
points instead of one -- and any lever held back in AMIP for fear of OVERCOOLING the SO
should be expected to overcool MORE coupled, not less.

SEASON MATTERS AND THE ANNUAL MEAN HIDES IT.  DMS and S4 agree on annual SO SW CRE by
coincidence and disagree by season: DJF -7.25 against -3.72.  DJF is where the deficit
lives (+10.86) because it carries the shortwave -- 460 W/m^2 incoming against 87 in JJA.
Scoring this pair on the annual mean alone would repeat an error the report already
withdrew, so DJF terms are scored explicitly.

AMIP ANCHORS (U1-T3, ONE year, so far weaker statistics than S4's 44-year anchor -- do
not read the ratios to more than one significant figure):
    SO 45-65S SW CRE   ANN -2.25   DJF -7.25

THRESHOLDS FIRST, from 11E's own interannual scatter over the same 20 years,
1.96*sd*sqrt(2/n).  Both arms branch from the one 1350 state, so this is a paired test.

A TRAP THIS AVOIDS.  IFS TOA fluxes are accumulated J/m^2 and must be divided by the
accumulation period; forgetting it makes every radiative number ~3600x too large.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
ACC = 3600.0

PAIR = [('11E base', f'{R092}/Tuning_test_11E_swemin15_K1'),
        ('11F +DMS', f'{R092}/Tuning_test_11F_dmsccn166')]

# 11F was cancelled at 20 years; this is everything it has.  Below the campaign's
# 30-year coupled minimum, which the threshold arithmetic handles honestly -- a
# shorter window simply raises the bar a difference has to clear.
CLEAN = list(range(1350, 1370))
DIRTY = []
TROPICS = (-20.0, 20.0)

SIB = (55.0, 75.0, 60.0, 180.0)
SO = (-65.0, -45.0)
DJF, JJA = [12, 1, 2], [6, 7, 8]
FLUX = ('tsr', 'ttr', 'tsrc', 'ttrc', 'ssr', 'str')


def years_avail(root, var):
    D = f'{root}/outdata/oifs'
    if not os.path.isdir(D):
        return []
    ys = []
    for f in os.listdir(D):
        if f.startswith(f'atm_remapped_1m_{var}_') and f.endswith('.nc'):
            try:
                ys.append(int(f.split('_')[-1].split('-')[0]))
            except ValueError:
                pass
    return sorted(ys)


def load_year(root, var, y):
    """(12, nlat, nlon) for one year, de-accumulated, plus lat/lon."""
    D = f'{root}/outdata/oifs'
    div = ACC if var in FLUX else 1.0
    for f in (f'{D}/atm_remapped_1m_{var}_{y}-{y}.nc',
              f'{D}/atm_remapped_1m_{var}_1m_{y}-{y}.nc'):
        if os.path.exists(f):
            with xr.open_dataset(f, decode_times=False) as d:
                a = d[var].values / div
                lat, lon = d['lat'].values, d['lon'].values
            return (a, lat, lon) if a.shape[0] == 12 else (None, None, None)
    return None, None, None


def gmean(f2d, lat):
    w = np.cos(np.deg2rad(lat))
    return float(np.average(f2d.mean(axis=1), weights=w))


def boxmean(f2d, lat, lon, box, lsm=None):
    la0, la1, lo0, lo1 = box
    ys = (lat >= la0) & (lat <= la1)
    xs = ((lon % 360) >= lo0) & ((lon % 360) <= lo1)
    sub = f2d[np.ix_(ys, xs)]
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape).copy()
    if lsm is not None:
        w = np.where(lsm[np.ix_(ys, xs)] >= 0.5, w, 0.0)
    return float(np.average(sub, weights=w)) if w.sum() else np.nan


def zband(f2d, lat, a, b):
    sel = (lat >= a) & (lat < b)
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(f2d[sel, :].mean(axis=1), weights=w))


def zband_ocean(f2d, lat, a, b, lsm):
    """Zonal band mean over OCEAN points only.

    The surface terms have to be ocean-masked: a band mean of SST or sea-ice
    cover that includes the Antarctic Peninsula and the southern tips of the
    continents is not a Southern Ocean number.
    """
    sel = (lat >= a) & (lat < b)
    w = np.broadcast_to(np.cos(np.deg2rad(lat[sel]))[:, None], f2d[sel, :].shape).copy()
    w = np.where(lsm[sel, :] < 0.5, w, 0.0)
    return float(np.average(f2d[sel, :], weights=w)) if w.sum() else np.nan


def metrics_for_year(root, y, lsm, lat, lon):
    """The scored quantities for a single year, or None if output is incomplete."""
    t2m, _, _ = load_year(root, '2t', y)
    if t2m is None:
        return None
    tsr, _, _ = load_year(root, 'tsr', y)
    ttr, _, _ = load_year(root, 'ttr', y)
    tsrc, _, _ = load_year(root, 'tsrc', y)
    stl2, _, _ = load_year(root, 'stl2', y)
    if tsr is None or ttr is None:
        return None
    d = [m - 1 for m in DJF]
    j = [m - 1 for m in JJA]
    out = {
        'net TOA [W/m2]': gmean((tsr + ttr).mean(axis=0), lat),
        'SO SW CRE [W/m2]': (zband((tsr - tsrc).mean(axis=0), lat, *SO)
                             if tsrc is not None else np.nan),
        'SO net TOA [W/m2]': zband((tsr + ttr).mean(axis=0), lat, *SO),
        # DJF carries the austral shortwave and is the season DMS actually works in;
        # the annual mean is where DMS and S4 agree by coincidence.
        'SO SW CRE DJF [W/m2]': (zband((tsr - tsrc)[d].mean(axis=0), lat, *SO)
                                 if tsrc is not None else np.nan),
        'tropics SW CRE [W/m2]': (zband((tsr - tsrc).mean(axis=0), lat, *TROPICS)
                                  if tsrc is not None else np.nan),
        'global T2m [C]': gmean(t2m.mean(axis=0), lat) - 273.15,
        'Siberia JJA T2m [C]': boxmean(t2m[j].mean(axis=0), lat, lon, SIB, lsm) - 273.15,
        'Siberia DJF T2m [C]': boxmean(t2m[d].mean(axis=0), lat, lon, SIB, lsm) - 273.15,
    }
    # The terms the coupled amplification actually lives in. Round 27 measured
    # 1.40x for DMS and found the coupled CRE response SMALLER than AMIP's --
    # the cloud adjusts away once the surface cools, and the amplification comes
    # from sea-ice albedo, which is CLEAR-SKY and invisible to CRE. So scoring an
    # SO lever on CRE alone measures the part that does not survive coupling.
    sst, _, _ = load_year(root, 'sst', y)
    ci, _, _ = load_year(root, 'ci', y)
    if tsrc is not None:
        out['SO clear-sky SW [W/m2]'] = zband(tsrc.mean(axis=0), lat, *SO)
    if sst is not None:
        out['SO SST [C]'] = zband_ocean(sst.mean(axis=0), lat, *SO, lsm) - 273.15
    if ci is not None:
        out['SO sea ice [%]'] = zband_ocean(ci.mean(axis=0), lat, *SO, lsm) * 100.0
    if stl2 is not None:
        out['Siberia DJF soil [C]'] = boxmean(stl2[d].mean(axis=0), lat, lon, SIB,
                                              lsm) - 273.15
        out['Siberia JJA soil [C]'] = boxmean(stl2[j].mean(axis=0), lat, lon, SIB,
                                              lsm) - 273.15
    return out


print(__doc__)
print('=' * 100)

avail = {}
for tag, root in PAIR:
    ys = years_avail(root, '2t')
    avail[tag] = ys
    print(f'  {tag:9s} {len(ys):3d} years on disk: {ys[0]}-{ys[-1]}' if ys
          else f'  {tag:9s} NO OUTPUT')
usable = [y for y in CLEAN if all(y in avail[t] for t, _ in PAIR)]
_excl = f', excluding {DIRTY[0]}-{DIRTY[-1]}' if DIRTY else ' (no exclusions needed)'
print(f'\n  matched clean window: {len(usable)} years '
      f'({usable[0]}-{usable[-1]}){_excl}')
if len(usable) < 30:
    print(f'  NOTE: {len(usable)} years is below the 30-year campaign minimum. The')
    print('        thresholds below already account for it -- a short window raises')
    print('        the bar rather than lowering it -- but an unresolved term here is')
    print('        NOT evidence of absence.')

series = {}
for tag, root in PAIR:
    lsm, lat, lon = load_year(root, 'lsm', usable[0])
    if lsm is not None and lsm.ndim == 3:
        lsm = lsm[0]
    rows = []
    for y in usable:
        m = metrics_for_year(root, y, lsm, lat, lon)
        if m:
            rows.append(m)
    series[tag] = rows
    print(f'  {tag:9s} scored {len(rows)} of {len(usable)} years')

keys = list(series[PAIR[0][0]][0].keys())

print('\n' + '=' * 100)
print(f'\nDETECTION THRESHOLDS, from 11E interannual scatter over the same '
      f'{len(usable)} years:  1.96*sd*sqrt(2/n)\n')
base_tag = PAIR[0][0]
thr = {}
for k in keys:
    v = np.array([r[k] for r in series[base_tag]])
    thr[k] = 1.96 * v.std(ddof=1) * np.sqrt(2.0 / len(v))
    print(f'  {k:24s} sd {v.std(ddof=1):8.4f}   threshold +-{thr[k]:.4f}')

print('\n' + '=' * 100)
print(f'\n11F minus 11E, {len(usable)} matched years.  * = resolved\n')
print(f'  {"metric":24s} {"11E":>10s} {"11F":>10s} {"diff":>11s} {"thr":>9s}')
verdict = {}
for k in keys:
    a = np.mean([r[k] for r in series[PAIR[0][0]]])
    b = np.mean([r[k] for r in series[PAIR[1][0]]])
    d = b - a
    sig = abs(d) > thr[k]
    verdict[k] = (a, b, d, sig)
    print(f'  {k:24s} {a:10.3f} {b:10.3f} {d:+10.3f}{"*" if sig else " "} '
          f'{thr[k]:9.3f}')

try:
    with xr.open_dataset(CERESF) as cds:
        clat = cds['lat'].values
        csel = (clat >= SO[0]) & (clat < SO[1])
        cw = np.cos(np.deg2rad(clat[csel]))
        so_obs = float(np.average(
            cds['toa_cre_sw_clim'].values.mean(axis=0)[csel, :].mean(axis=1), weights=cw))
    a, b, _, _ = verdict['SO SW CRE [W/m2]']
    print(f'\n  CERES SO SW CRE {so_obs:.2f}:  11E is {a - so_obs:+.2f} from it, '
          f'11F {b - so_obs:+.2f}')
except Exception as exc:
    print(f'\n  (CERES anchor unavailable: {exc})')

print('\n' + '=' * 100)
print('\nPRE-REGISTERED READING\n')
so_ann = verdict['SO SW CRE [W/m2]']
so_djf = verdict.get('SO SW CRE DJF [W/m2]')
trop = verdict.get('tropics SW CRE [W/m2]')
AMIP_ANN, AMIP_DJF = -2.25, -7.25          # U1-T3, 1 yr, SO 45-65S

print(f'  SO SW CRE ANN  {so_ann[0]:+.3f} -> {so_ann[1]:+.3f}  ({so_ann[2]:+.3f}, '
      f'{"resolved" if so_ann[3] else "NOT resolved"});  AMIP {AMIP_ANN:+.2f}')
if so_djf:
    print(f'  SO SW CRE DJF  {so_djf[0]:+.3f} -> {so_djf[1]:+.3f}  ({so_djf[2]:+.3f}, '
          f'{"resolved" if so_djf[3] else "NOT resolved"});  AMIP {AMIP_DJF:+.2f}')
if trop:
    print(f'  tropics SW CRE {trop[2]:+.3f} ({"resolved" if trop[3] else "ns"})'
          '   <- the cost that disqualifies DMS as a knob, regardless of this result')
print()
print('  COUPLED / AMIP RATIO -- the actual question:')
for lbl, v, amip in (('ANN', so_ann, AMIP_ANN), ('DJF', so_djf, AMIP_DJF)):
    if v:
        print(f'    {lbl}  {v[2] / amip:.2f}x      (S4 gave 0.28-0.31x; DMS was '
              f'predicted 1.40x)')
print()
for k in ('SO clear-sky SW [W/m2]', 'SO SST [C]', 'SO sea ice [%]'):
    if k in verdict:
        a, b, dd, sg = verdict[k]
        print(f'  {k:24s} {dd:+8.3f} ({"resolved" if sg else "ns"})')
print()
ice = verdict.get('SO sea ice [%]')
sst = verdict.get('SO SST [C]')
print('  THE SCREENING RULE THIS TESTS:')
if ice and sst:
    print(f'    SO SST {sst[2]:+.3f} K, sea ice {ice[2]:+.3f} pp.  S4 managed -0.264 K '
          'and +0.595 pp')
    print('    and DAMPED to 0.28x. If this arm moves the surface much harder AND')
    print('    amplifies, the discriminator is confirmed: a lever feared for overcooling')
    print('    the SO in AMIP will overcool MORE coupled, not less.')
print()
for k in ('Siberia JJA T2m [C]', 'Siberia DJF T2m [C]'):
    if k in verdict:
        a, b, dd, sg = verdict[k]
        print(f'  {k:24s} {dd:+8.3f} ({"resolved" if sg else "ns"})')
