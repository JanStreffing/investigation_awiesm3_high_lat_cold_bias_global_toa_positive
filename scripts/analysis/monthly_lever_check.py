"""Month-by-month response of a lever, for testing melt-timing predictions.

Round 13 predicted that raising the snow-cover critical depth (RQSNCR 1/10 -> 1/30)
would recover June's albedo term and do LITTLE in July and August, because by then
the surface is snow-free and the residual SW deficit is cloud. That prediction is
only testable per month -- a JJA mean averages the claim away.

Reports, per calendar month, on the Siberian land box:
  * T2m change vs control
  * surface net SW change vs control
  * surface albedo change (derived as 1 - ssr/ssrd, so it is the ALL-SKY surface
    albedo actually seen by the radiation, not a tile-weighted diagnostic)

with a per-month noise floor from the same run x year ANOVA used elsewhere, because
a single month is noisier than a season and monthly deltas are easy to over-read.

Usage:  python3 monthly_lever_check.py [label ...]
        defaults to the H-series and its parents.
"""
import numpy as np, xarray as xr, os, sys, warnings
warnings.filterwarnings('ignore')

from runs import RUNS, RT, LSMF, Y0, Y1

YEARS = list(range(Y0, Y1 + 1))
BOX = ((55, 75), (60, 180))
ACC = 3600.0
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

WANT = sys.argv[1:] or ['F4 rsmin1000', 'G1 F4+D2b', 'H1 snowcr30', 'H2 G1+snowcr']

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values
yi = None


def box_mean(a2d, lat, lon):
    global yi
    ys = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
    l180 = ((lon + 180) % 360) - 180
    xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    sub = a2d[np.ix_(np.where(ys)[0], np.where(xs)[0])]
    L = lsm[np.ix_(np.where(ys)[0], np.where(xs)[0])]
    m = L > 0.5
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape)
    return np.average(sub[m], weights=w[m])


def field(run, var, y):
    f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
    if not os.path.exists(f):
        return None, None, None
    ds = xr.open_dataset(f)
    a = ds[var].values
    lat, lon = ds[var].lat.values, ds[var].lon.values
    ds.close()
    return a, lat, lon


def per_year_monthly(run):
    """(t2m, ssr, albedo) each [n_year, 12] on the Siberian box; None if incomplete."""
    T, S, A = [], [], []
    for y in YEARS:
        t, lat, lon = field(run, '2t', y)
        sr, _, _ = field(run, 'ssr', y)
        sd, _, _ = field(run, 'ssrd', y)
        if t is None or sr is None:
            return None
        T.append([box_mean(t[m], lat, lon) - 273.15 for m in range(12)])
        S.append([box_mean(sr[m] / ACC, lat, lon) for m in range(12)])
        if sd is not None:
            # all-sky surface albedo = 1 - net/down, computed on the box means so
            # it is the energetically-weighted albedo, not an area-mean of ratios
            A.append([1.0 - box_mean(sr[m] / ACC, lat, lon) /
                      max(box_mean(sd[m] / ACC, lat, lon), 1e-6) for m in range(12)])
    return np.array(T), np.array(S), (np.array(A) if A else None)


data, labs = {}, []
for lab, run in RUNS:
    if lab != 'control' and lab not in WANT:
        continue
    r = per_year_monthly(run)
    if r is None:
        print(f'  !! {lab}: incomplete, skipped'); continue
    labs.append(lab); data[lab] = r
if 'control' not in labs:
    sys.exit('control missing')


def report(idx, name, unit, scale=1.0):
    print(f'\n=== {name} vs control, Siberia land [{unit}] ===')
    X = np.array([data[l][idx] for l in labs])          # [run, year, month]
    ctl = labs.index('control')
    hdr = '  ' + ' '.join(f'{m:>8s}' for m in MON)
    print(f"  {'':16s}" + hdr)
    thr = []
    for m in range(12):
        Y = X[:, :, m]
        mu = Y.mean(); a = Y.mean(1) - mu; g = Y.mean(0) - mu
        eps = Y - (mu + a[:, None] + g[None, :])
        sd = np.sqrt((eps ** 2).sum() / ((Y.shape[0] - 1) * (Y.shape[1] - 1)))
        thr.append(1.96 * sd * np.sqrt(2.0 / Y.shape[1]))
    print(f"  {'95% threshold':16s}  " + ' '.join(f'{t*scale:>8.3f}' for t in thr))
    print(f"  {'control':16s}  " + ' '.join(f'{X[ctl,:,m].mean()*scale:>8.2f}' for m in range(12)))
    for i, l in enumerate(labs):
        if l == 'control':
            continue
        d = [X[i, :, m].mean() - X[ctl, :, m].mean() for m in range(12)]
        print(f'  {l:16s}  ' + ' '.join(
            f'{v*scale:>+7.3f}{"*" if abs(v) > thr[m] else " "}' for m, v in enumerate(d)))


report(0, 'T2m', 'K')
report(1, 'surface net SW', 'W/m2')
if data[labs[0]][2] is not None:
    report(2, 'all-sky surface albedo', 'x100', scale=100.0)
print('\n  * = clears that month\'s own 95% detection threshold.')
