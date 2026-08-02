"""Siberian T2m by season for every run, with a per-season noise floor.

The campaign evaluated JJA only, for eleven rounds. That is how B5's -0.72 K
winter cooling survived unnoticed until round 12 -- and the coupled model's
original complaint is a *cold-season* bias, so a lever that fixes summer by
breaking winter is worse than useless.

This script runs the same run x year ANOVA as noise_floor.py, but separately
for each season, so every seasonal delta comes with its own detection
threshold. The noise floor is NOT the same in each season: Siberian winter
interannual variability is several times summer's, so a winter delta needs to
be much larger before it means anything.

Same box, land mask, window and control as eval_round10_A.py, so the JJA
column here must reproduce that table's JJA row exactly. It is a deliberate
cross-check -- if the two disagree, one of them is wrong.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

from runs import RUNS, RT, LSMF, Y0, Y1

YEARS = list(range(Y0, Y1 + 1))
BOX = ((55, 75), (60, 180))          # Siberia land, as everywhere else

# 0-based month indices, matching eval_round10_A.py's convention (JJA = 5,6,7).
# DJF uses Dec of the SAME calendar year rather than crossing the year
# boundary: with 44 years the difference is a 1/44 edge effect on one month,
# far below the winter noise floor, and it keeps every season on one file.
SEASONS = {'DJF': [11, 0, 1], 'MAM': [2, 3, 4], 'JJA': [5, 6, 7], 'SON': [8, 9, 10]}

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def box_mean(a2d, lat, lon):
    yi = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
    l180 = ((lon + 180) % 360) - 180
    xi = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    sub = a2d[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    L = lsm[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    m = L > 0.5
    w = np.broadcast_to(np.cos(np.deg2rad(lat[yi]))[:, None], sub.shape)
    return np.average(sub[m], weights=w[m])


def per_year_seasons(run):
    """{season: array over years} of Siberian mean T2m; None if any year missing."""
    out = {s: [] for s in SEASONS}
    for y in YEARS:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_2t_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            return None
        ds = xr.open_dataset(f)
        a = ds['2t'].values
        lat, lon = ds['2t'].lat.values, ds['2t'].lon.values
        ds.close()
        mon = [box_mean(a[m], lat, lon) for m in range(12)]
        for s, idx in SEASONS.items():
            out[s].append(np.mean([mon[i] for i in idx]))
    return {s: np.array(v) - 273.15 for s, v in out.items()}


labs, data = [], {s: [] for s in SEASONS}
for lab, run in RUNS:
    r = per_year_seasons(run)
    if r is None:
        print(f'  !! {lab}: incomplete, skipped'); continue
    labs.append(lab)
    for s in SEASONS:
        data[s].append(r[s])

ctl = labs.index('control')
thr, delta = {}, {}
print(f'\nSiberia 55-75N 60-180E land, T2m by season. {len(labs)} runs x '
      f'{len(YEARS)} yr ({Y0}-{Y1}). Run x year ANOVA per season.\n')

for s in SEASONS:
    X = np.array(data[s]); nr, ny = X.shape
    mu = X.mean(); a = X.mean(1) - mu; g = X.mean(0) - mu
    eps = X - (mu + a[:, None] + g[None, :])
    sd = np.sqrt((eps ** 2).sum() / ((nr - 1) * (ny - 1)))
    thr[s] = 1.96 * sd * np.sqrt(2.0 / ny)
    delta[s] = X.mean(1) - X[ctl].mean()
    print(f'  {s}: control {X[ctl].mean():6.2f} C   sd(eps)={sd:5.3f} K   '
          f'95% threshold +-{thr[s]:.3f} K')

print(f"\n  {'run':16s} " + ' '.join(f'{s:>17s}' for s in SEASONS))
for i, lab in enumerate(labs):
    row = ''
    for s in SEASONS:
        d = delta[s][i]
        mark = '*' if abs(d) > thr[s] and i != ctl else ' '
        row += f'  {d:+8.3f}{mark} {"":6s}' if False else f'   {d:+8.3f}{mark}      '
    print(f'  {lab:16s} {row}')
print('\n  * = clears that season\'s own 95% detection threshold.')

# The specific failure mode this script exists to catch: summer gain, winter loss.
print('\n  Levers that WARM JJA significantly but COOL DJF significantly:')
bad = [labs[i] for i in range(len(labs)) if i != ctl
       and delta['JJA'][i] > thr['JJA'] and delta['DJF'][i] < -thr['DJF']]
print('    ' + (', '.join(bad) if bad else 'none'))
