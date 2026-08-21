"""Does Raupach (1994) canopy roughness move the NH-WINTER land cold bias?

WHY.  land_bias_season_band.py localised the coupled land cold bias to winter: the
land-minus-ocean T2m gap is -1.37 K in DJF against -0.09 K in JJA, i.e. the land-specific
error is a cold-season one, consistent with the campaign's measured "+1.4 K too strong"
low-level inversion.  Round 16 (snow skin conductivity, ECE_LAMSK_SN 7 -> 15/25) was the
first candidate and is now scored NULL -- see lamsk_j_series_band.py; its pre-registered
falsifier said that if DJF does not respond the winter bias lies in boundary-layer mixing
instead.  Roughness length is the other term in that exchange.

THE LEVER.  As released, EC-Earth/OpenIFS derives z0 for high vegetation from a per-
vegetation-type lookup driven by the high-vegetation fraction, so z0 RISES with high-veg
LAI.  Raupach (1994) predicts the opposite: z0 FALLS as the canopy closes and smooths
with increasing LAI.  In a boreal winter -- low LAI -- the two schemes therefore differ in
the direction that matters, Raupach giving the rougher surface, more mechanical mixing,
and a weaker inversion.  Commit 1385d0c in lpj_guess_50yr_raupach computes z0 from FPC,
canopy height and a per-PFT frontal-area index and ships it as GUE_Z0HV -> Z0HVeg,
replacing the lookup, gated on the run-time switch ifraupachz0.

THE PAIR.  11J = 11I + ifraupachz0, nothing else, and both branch from the same 1350
state, so the difference is attributable to the roughness physics rather than to the
build (one binary, run-time gate).  Verified before scoring: 11J's guess.ins carries
"ifraupachz0 1" and fort.4 "ECE_CPL_LPJG_Z0 = .true.", 11I has neither.

Paired statistics throughout -- branched coupled arms, so the threshold is
1.96 * sd(annual differences) / sqrt(n), not the two-sample one.

READ THE SEASONS AND READ THE LAND-MINUS-OCEAN GAP.  A lever that warms land and ocean
alike has not touched the land bias; it has moved the global mean, which the TOA budget
already spends.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CTL = ('11I v2soil', f'{R}/Tuning_test_11I_v2soil')
ARM = ('11J +raupach', f'{R}/Tuning_test_11J_v2soil_raupach')
# Two windows: the campaign's standard decade, plus the earlier one as a check that the
# signal is not a single-decade fluctuation of a drifting coupled pair.
WINDOWS = [(1380, 1389), (1360, 1379)]
SEASONS = [('DJF', [12, 1, 2]), ('MAM', [3, 4, 5]),
           ('JJA', [6, 7, 8]), ('SON', [9, 10, 11]), ('ANN', list(range(1, 13)))]
BANDS = [('60-90N', 60, 90), ('30-60N', 30, 60), ('0-30N', 0, 30),
         ('30S-0', -30, 0), ('60-30S', -60, -30)]


def grid(path):
    f = sorted(glob.glob(f'{path}/outdata/oifs/atm_remapped_1m_lsm_*.nc'))[0]
    with xr.open_dataset(f, decode_times=False) as d:
        m = np.squeeze(d['lsm'].values)
        if m.ndim == 3:
            m = m[0]
        la = np.squeeze(d['lat'].values)
    return m, la


def load(path, y0, y1):
    out = []
    for y in range(y0, y1 + 1):
        f = f'{path}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc'
        if not os.path.exists(f):
            return None
        with xr.open_dataset(f, decode_times=False) as d:
            a = np.squeeze(d['2t'].values)
        if a.shape[0] != 12:
            raise SystemExit(f'{f}: {a.shape[0]} steps, expected 12')
        out.append(a)
    a = np.stack(out)
    return a - 273.15 if np.nanmean(a) > 100 else a


msk, lat = grid(CTL[1])
land, ocean = msk > 0.5, msk <= 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], msk.shape).copy()


def am(field, sel):
    m = sel & np.isfinite(field)
    return float(np.average(field[m], weights=W[m])) if m.any() else np.nan


print(f'Raupach z0: {ARM[0]} minus {CTL[0]} (coupled, paired). '
      f'* = |delta| > paired 95% threshold.')

for y0, y1 in WINDOWS:
    b = load(CTL[1], y0, y1)
    a = load(ARM[1], y0, y1)
    if a is None or b is None:
        print(f'\n[{y0}-{y1}] incomplete, skipped')
        continue
    print(f'\n=== {y0}-{y1} ({y1 - y0 + 1} yr) ===')
    print(f'{"band":>8} | ' + ' | '.join(f'{s:>15}' for s, _ in SEASONS))
    rows = BANDS + [('ALL LAND', None, None), ('OCEAN', None, None), ('GAP', None, None)]
    for bname, lo, hi in rows:
        cells = []
        for _, mons in SEASONS:
            mi = [m - 1 for m in mons]
            dif = [a[y][mi].mean(0) - b[y][mi].mean(0) for y in range(a.shape[0])]
            if bname == 'OCEAN':
                dy = np.array([am(d, ocean) for d in dif])
            elif bname == 'ALL LAND':
                dy = np.array([am(d, land) for d in dif])
            elif bname == 'GAP':
                dy = np.array([am(d, land) - am(d, ocean) for d in dif])
            else:
                sel = land & np.broadcast_to((lat >= lo) & (lat < hi), msk.shape[::-1]).T \
                    if False else land & np.broadcast_to(
                        ((lat >= lo) & (lat < hi))[:, None], msk.shape)
                dy = np.array([am(d, sel) for d in dif])
            m_ = float(np.nanmean(dy))
            thr = 1.96 * float(np.nanstd(dy, ddof=1)) / np.sqrt(np.isfinite(dy).sum())
            cells.append(f'{m_:+7.3f} {"*" if abs(m_) > thr else " "}({thr:.2f})')
        pre = '' if bname not in ('GAP',) else ''
        print(f'{bname:>8} | ' + ' | '.join(f'{c:>15}' for c in cells))
