"""Round 16, scored at last: is snow skin conductivity the NH-winter LAND lever?

WHY.  land_bias_season_band.py resolved the coupled land cold bias to a WINTER,
NH-extratropical, land-specific error: the land-minus-ocean gap is -1.37 K in DJF and
only -0.09 K in JJA.  Round 16 (2026-08-04) raised the exposed-snow skin conductivity
ZSNOW from 7 to 15 (J1) and 25 (J2) on top of I1 and predicted "DJF +1...+3 K, JJA
barely moves" -- exactly that shape.  Both runs completed 48 years and were never
scored; the round is still marked "in flight" in RUNS_AND_PARAMETERS.md.

The staged artefact was checked before writing this: libsurf.SP.so (NOT libarpifs) in
both run directories carries ECE_LAMSK_SN, built 2026-08-04, and fort.4 has 15 / 25.
So this is not another silent no-op like W1.

PAIRING.  These are AMIP arms under identical prescribed SST, started from the same
state, so the year-by-year difference J-I1 is paired and the threshold is the paired
one: 1.96 * sd(annual differences) / sqrt(n).  Using the unpaired two-sample threshold
here would be far too lax and would promote noise.

READ THE SEASONS, NOT THE ANNUAL MEAN, and read the land-minus-ocean gap rather than
land alone -- a lever that warms land and ocean equally has not touched the land bias,
it has just shifted the global mean, which is already spoken for by the TOA budget.
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

from runs import RT, LSMF, Y0, Y1

BASE = ('I1', 'amip_I1_scf')
ARMS = [('J1 lamsk15', 'amip_J1_lamsk15'), ('J2 lamsk25', 'amip_J2_lamsk25')]
SEASONS = [('DJF', [12, 1, 2]), ('MAM', [3, 4, 5]),
           ('JJA', [6, 7, 8]), ('SON', [9, 10, 11]), ('ANN', list(range(1, 13)))]
BANDS = [('60-90N', 60, 90), ('30-60N', 30, 60), ('0-30N', 0, 30),
         ('30S-0', -30, 0), ('60-30S', -60, -30)]


def lsm():
    with xr.open_dataset(LSMF, decode_times=False) as d:
        v = [k for k in d.data_vars if k.lower() in ('lsm', 'var172')][0]
        m = np.squeeze(d[v].values)
        if m.ndim == 3:      # lsm is written on the monthly axis; it is constant
            m = m[0]
        la = np.squeeze(d[[c for c in d[v].dims if 'lat' in c][0]].values)
    return m, la


def load(run):
    """(year, month, ny, nx) T2m in C for the shared window."""
    out = []
    for y in range(Y0, Y1 + 1):
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_2t_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            return None
        with xr.open_dataset(f, decode_times=False) as d:
            v = [k for k in d.data_vars if k.lower() in ('2t', 't2m', 'var167')][0]
            a = np.squeeze(d[v].values)
        if a.shape[0] != 12:
            raise SystemExit(f'{f}: {a.shape[0]} steps, expected 12')
        out.append(a)
    a = np.stack(out)
    return a - 273.15 if np.nanmean(a) > 100 else a


def area_mean(field, mask):
    """Equal-area: the remapped grid is regular lat-lon, so cos-lat weight."""
    m = mask & np.isfinite(field)
    if not m.any():
        return np.nan
    return float(np.average(field[m], weights=W[m]))


msk, lat = lsm()
land = msk > 0.5
ocean = msk <= 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], msk.shape).copy()

base = load(BASE[1])
if base is None:
    raise SystemExit(f'base {BASE[1]} incomplete over {Y0}-{Y1}')

print(f'Round 16 scored: ECE_LAMSK_SN (exposed-snow skin conductivity, default 7)')
print(f'AMIP, paired against {BASE[0]}, {Y0}-{Y1} ({Y1 - Y0 + 1} yr). '
      f'* = |delta| > paired 95% threshold.\n')

for label, run in ARMS:
    arm = load(run)
    if arm is None:
        print(f'{label}: incomplete, skipped\n')
        continue
    print(f'=== {label} minus {BASE[0]} ===')
    print(f'{"band":>8} | ' + ' | '.join(f'{s:>16}' for s, _ in SEASONS))
    for bname, lo, hi in BANDS + [('ALL LAND', -90, 90), ('OCEAN', -90, 90)]:
        if bname == 'OCEAN':
            sel = ocean
        else:
            band = (lat >= lo) & (lat < hi)
            sel = land & np.broadcast_to(band[:, None], msk.shape)
        cells = []
        for _, mons in SEASONS:
            mi = [m - 1 for m in mons]
            # DJF: keep it simple and consistent with the rest of the campaign --
            # calendar-year DJF, same convention in both arms, so the pairing holds.
            dy = np.array([area_mean(arm[y][mi].mean(0) - base[y][mi].mean(0), sel)
                           for y in range(arm.shape[0])])
            d = float(np.nanmean(dy))
            thr = 1.96 * float(np.nanstd(dy, ddof=1)) / np.sqrt(np.isfinite(dy).sum())
            cells.append(f'{d:+7.3f} {"*" if abs(d) > thr else " "}({thr:.2f})')
        print(f'{bname:>8} | ' + ' | '.join(f'{c:>16}' for c in cells))
    # the number that actually matters
    print(f'{"":>8} | land-minus-ocean gap, by season:')
    gaps = []
    for _, mons in SEASONS:
        mi = [m - 1 for m in mons]
        dy = np.array([area_mean(arm[y][mi].mean(0) - base[y][mi].mean(0), land)
                       - area_mean(arm[y][mi].mean(0) - base[y][mi].mean(0), ocean)
                       for y in range(arm.shape[0])])
        d = float(np.nanmean(dy))
        thr = 1.96 * float(np.nanstd(dy, ddof=1)) / np.sqrt(np.isfinite(dy).sum())
        gaps.append(f'{d:+7.3f} {"*" if abs(d) > thr else " "}({thr:.2f})')
    print(f'{"GAP":>8} | ' + ' | '.join(f'{g:>16}' for g in gaps))
    print()
