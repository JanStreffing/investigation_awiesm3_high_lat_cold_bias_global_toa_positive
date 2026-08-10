"""How much does the coupled Southern Ocean amplify an AMIP forcing?

THE QUESTION.  Every candidate configuration has been judged against the AMIP Southern
Ocean gap -- 7.85 W/m2 of SW CRE against CERES, 10.86 in DJF -- and on that scale they all
look inadequate: the best reaches about a quarter.  But AMIP holds SST and sea ice fixed,
which suppresses the feedback that makes a Southern Ocean cloud change matter.  Coupled,
more reflection cools SST, ice expands, albedo rises, and that cools further.  So the AMIP
number is a LOWER BOUND on the coupled response, and the amount a lever must deliver in
AMIP is correspondingly less than the full gap.

THE MEASUREMENT, and it costs nothing.  11E and 11F are a clean coupled pair over their
common 20 years (1350-1369): identical in every respect except ECE_DMS_CCN_SENS, 0 against
166.  That is a KNOWN Southern Ocean perturbation -- the same lever measured in AMIP as
U1 minus T3, giving SO SW CRE -2.250 and SO net -1.976 with SST prescribed.  Comparing the
two gives the amplification factor directly:

    amplification  =  (coupled d SO net TOA)  /  (AMIP d SO net TOA)

and the state response -- how much SST and sea ice actually moved -- says whether the
feedback is engaging at all.

WHY 11F EXISTS TO BE USED.  It was cancelled at 20 of 50 years on 2026-08-10 as too
expensive for what it does, after its lever was superseded.  The 20 years are clean
(both legs CNT0) and this is the one question they answer that nothing else can.

THRESHOLDS ARE COMPUTED FIRST, from 11E's own interannual scatter over the same 20 years,
as 1.96*sd*sqrt(2/20) for a 20-year pair.  Coupled fields are far noisier than AMIP ones
because the ocean has its own variability, and the campaign has already been burned by
quoting a coupled delta without one: the measured detection threshold on a 10-year coupled
pair is +-0.75 K for soil and +-1.25 K for DJF T2m.

CAVEAT.  Twenty years is short for sea ice, which adjusts on decadal timescales, so a null
on ice area is a bound rather than evidence of no feedback.  The TOA and SST responses are
faster and should be readable.
"""
import os, glob
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

B = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CTL, ARM = 'Tuning_test_11E_swemin15_K1', 'Tuning_test_11F_dmsccn166'
Y0, Y1 = 1350, 1369
ACC = 3600.0            # coupled OpenIFS radiative accumulation, same as AMIP
SO = (-65.0, -45.0)
# the same lever measured in AMIP, U1 - T3, year 1870
AMIP = {'SO SW CRE': -2.250, 'SO LW CRE': 0.389, 'SO net TOA': -1.976}

print(__doc__)
print('=' * 100)


def oifs_years(run, var):
    """Per-year SO-band mean of an OpenIFS monthly field."""
    out = []
    for y in range(Y0, Y1 + 1):
        f = f'{B}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            f2 = f'{B}/{run}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
            if not os.path.exists(f2):
                return None, None
            f = f2
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values
            lat = d['lat'].values
        out.append(a)
    return np.array(out), lat


def band(a3, lat, lo=SO[0], hi=SO[1]):
    s = (lat >= lo) & (lat < hi)
    w = np.cos(np.deg2rad(lat[s]))
    return float(np.average(a3.mean(axis=0)[s, :].mean(axis=1), weights=w))


def series(run):
    """Per-year SO values of the radiative and surface metrics."""
    got = {}
    for v in ('tsr', 'tsrc', 'ttr', 'ttrc', 'sst', 'ci'):
        a, lat = oifs_years(run, v)
        if a is None:
            return None
        got[v] = a / (ACC if v in ('tsr', 'tsrc', 'ttr', 'ttrc') else 1.0)
        got['lat'] = lat
    lat = got['lat']
    n = got['tsr'].shape[0]
    rows = {k: [] for k in ('SO SW CRE', 'SO LW CRE', 'SO net TOA', 'SO SST', 'SO ice frac')}
    for i in range(n):
        sw = got['tsr'][i] - got['tsrc'][i]
        lw = got['ttr'][i] - got['ttrc'][i]
        net = got['tsr'][i] + got['ttr'][i]
        rows['SO SW CRE'].append(band(sw, lat))
        rows['SO LW CRE'].append(band(lw, lat))
        rows['SO net TOA'].append(band(net, lat))
        rows['SO SST'].append(band(got['sst'][i], lat))
        rows['SO ice frac'].append(band(got['ci'][i], lat))
    return {k: np.array(v) for k, v in rows.items()}


c, a = series(CTL), series(ARM)
if c is None or a is None:
    raise SystemExit('one of the runs lacks the full 1350-1369 output')

n = len(c['SO SW CRE'])
print(f'11E (control) and 11F (DMS S=166), {n} common years {Y0}-{Y1}\n')
print(f'  {"metric":14s} {"11E":>10s} {"11F":>10s} {"delta":>9s} {"thr(20yr)":>10s} '
      f'{"AMIP delta":>11s} {"amplification":>14s}')
print('  ' + '-' * 84)
for k in ('SO SW CRE', 'SO LW CRE', 'SO net TOA', 'SO SST', 'SO ice frac'):
    d = a[k].mean() - c[k].mean()
    thr = 1.96 * c[k].std(ddof=1) * np.sqrt(2.0 / n)
    am = AMIP.get(k)
    amp = f'{d/am:14.2f}' if am and abs(am) > 1e-9 else ' ' * 14
    amps = f'{am:11.3f}' if am else ' ' * 11
    star = '*' if abs(d) > thr else ' '
    print(f'  {k:14s} {c[k].mean():10.3f} {a[k].mean():10.3f} {d:9.3f}{star} {thr:10.3f} '
          f'{amps} {amp}')

print("\n  '*' = above the 20-year pair threshold from 11E's own interannual scatter.")
dnet = a['SO net TOA'].mean() - c['SO net TOA'].mean()
thrnet = 1.96 * c['SO net TOA'].std(ddof=1) * np.sqrt(2.0 / n)
print()
if abs(dnet) <= thrnet:
    print('  The coupled net TOA response is NOT resolved at 20 years, so the')
    print('  amplification factor is a bound, not a measurement.  Note this does not')
    print('  mean the feedback is absent -- coupled TOA is much noisier than AMIP TOA,')
    print('  and the SST/ice rows above are the more direct evidence.')
else:
    f = dnet / AMIP['SO net TOA']
    print(f'  Coupled net TOA response is {f:.2f}x the AMIP one.')
    if f > 1.2:
        print('  => the feedback AMPLIFIES, so the AMIP gap overstates what a lever must')
        print('     deliver, and candidates should not be judged against the full 7.85.')
    elif f < 0.8:
        print('  => the coupled response is WEAKER than AMIP, which would mean the ocean')
        print('     is absorbing the perturbation rather than amplifying it.')
    else:
        print('  => roughly one-to-one; no strong amplification at this length.')

print('\n  Sign guide: SO SST and ice fraction are the state response.  Cooling SST and')
print('  EXPANDING ice (positive d ice) is the feedback engaging in the intended')
print('  direction; the SO bias in coupled mode is a WARM one with too little ice.')
