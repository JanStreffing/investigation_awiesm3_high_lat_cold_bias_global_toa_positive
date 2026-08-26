"""Task 2: what meridional shape did the adopted stack give the Southern Ocean?

WHY.  Every SH-wide brightening lever -- DMS, RCL_KK_CLOUD_NUM_SEA -- works at the pole
and dies on 60-45S, which sits at -3.70 W/m2 DJF, overcorrected.  We put it there.  If one
adopted lever is EQUATORWARD-weighted, contributing more to that overcorrection than to
the polar fix, backing it off frees headroom for a lever family that is otherwise usable.

So decompose the stack step by step and ask, for each: what is the ratio of its effect at
90-60S ocean to its effect at 60-45S ocean?  A ratio above 1 is poleward-weighted and
worth more; below 1 is equatorward-weighted and is spending the budget in the wrong place.

WHY 30 YEARS AND PAIRED.  A first attempt differenced the 10-year arm table and got
CONTRADICTORY signs for the same lever: RSBLB read +1.57 at 90-60S from 11Q-11N and -3.55
from 11R-11P.  These are branched coupled arms, so the difference is paired and the
threshold is 1.96*sd/sqrt(n); at n=10 on this band that threshold is around 1 W/m2, which
is why the two disagreed.  Thirty years is the minimum that separates them.

PAIRS, each one namelist change apart and branched from the same 1350 state:
    S4  / RCL_INPPMIN   11G - 11E
    LX4 / RSNOWLIN2     11N - 11G
    RSBLB               11Q - 11N   (1850)  and  11R - 11P  (1990)
Both RSBLB pairs are scored because the 1850/1990 duplicate is the only internal check
available on a lever whose two estimates previously disagreed in sign.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
Y0, Y1 = 1360, 1389
PAIRS = [('S4 / INPPMIN', f'{R}/Tuning_test_11G_inppmin50k', f'{R}/Tuning_test_11E_swemin15_K1'),
         ('LX4 / RSNOWLIN2', f'{R}/11N', f'{R}/Tuning_test_11G_inppmin50k'),
         ('RSBLB (1850)', f'{R}/11Q', f'{R}/11N'),
         ('RSBLB (1990)', f'{R}/11R', f'{R}/11P')]
DJF, ANN = [11, 0, 1], list(range(12))

f = sorted(glob.glob(f'{R}/11P/outdata/oifs/atm_remapped_1m_lsm_*.nc'))[0]
with xr.open_dataset(f, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values)
    m = m[0] if m.ndim == 3 else m
    lat = np.squeeze(d['lat'].values)
ocean = m <= 0.5
RE = 6.371e6
A = np.broadcast_to((RE**2 * np.cos(np.deg2rad(lat)) * np.deg2rad(abs(lat[1] - lat[0]))
                     * 2 * np.pi / m.shape[1])[:, None], m.shape)
POLAR = ocean & np.broadcast_to((lat < -60)[:, None], m.shape)
SUBANT = ocean & np.broadcast_to(((lat >= -60) & (lat < -45))[:, None], m.shape)


def yrs(path, var):
    out = []
    for y in range(Y0, Y1 + 1):
        p = f'{path}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(p):
            return None
        with xr.open_dataset(p, decode_times=False) as d:
            out.append(np.squeeze(d[[k for k in d.data_vars
                                     if 'bnds' not in k and 'bounds' not in k][0]].values))
    return np.stack(out)


def am(f_, s):
    k = s & np.isfinite(f_)
    return float(np.average(f_[k], weights=A[k]))


def stat(dy):
    dy = np.asarray(dy)
    return float(np.nanmean(dy)), 1.96 * float(np.nanstd(dy, ddof=1)) / np.sqrt(len(dy))


print(__doc__)
print('=' * 96)
print(f'SW CRE change per lever, OCEAN-masked, paired, {Y0}-{Y1} ({Y1-Y0+1} yr) [W/m2]\n')
print(f'{"lever":>17} | {"polar 90-60S DJF":>20} {"subAnt 60-45S DJF":>20} | {"ratio":>7}')
print('-' * 96)
for lab, ap, cp in PAIRS:
    a1, a2 = yrs(ap, 'tsr'), yrs(ap, 'tsrc')
    c1, c2 = yrs(cp, 'tsr'), yrs(cp, 'tsrc')
    if any(x is None for x in (a1, a2, c1, c2)):
        print(f'{lab:>17} | incomplete'); continue
    n = min(a1.shape[0], c1.shape[0])
    cells = []
    for sel in (POLAR, SUBANT):
        mu, th = stat([am(((a1[i][DJF].mean(0) - a2[i][DJF].mean(0))
                          - (c1[i][DJF].mean(0) - c2[i][DJF].mean(0))) / 3600., sel)
                       for i in range(n)])
        cells.append((mu, th))
    r = cells[0][0] / cells[1][0] if abs(cells[1][0]) > 1e-6 else np.nan
    print(f'{lab:>17} | ' + ' '.join(f'{mu:+14.2f}{"*" if abs(mu) > th else " "}({th:4.2f})'
                                     for mu, th in cells) + f' | {r:7.2f}')
print('\nratio > 1 = POLEWARD-weighted (buys the target cheaply)')
print('ratio < 1 = EQUATORWARD-weighted (spends the 60-45S budget for little polar gain)')
