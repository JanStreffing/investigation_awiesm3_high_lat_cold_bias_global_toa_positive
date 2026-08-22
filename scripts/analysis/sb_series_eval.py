"""SB series, preliminary: does lowering RSBLB weaken the winter inversion?

WHY THIS FORM.  sb_run_length.py showed the metric to score is the INVERSION
(T925 - T2m), not the screen temperature: the paired annual sd at 60-90N is 0.611 K on
the inversion against 1.377 K on T2m, because differencing two levels in one column
cancels the synoptic variance.  That is a factor ~5 in years-to-decide, which is what
makes a 16-year read worth taking at all.

PAIRING.  SB0/SB1/SB2 are AMIP arms under identical prescribed SST from the same start,
differing in one namelist number, so the year-by-year difference is paired and the
threshold is 1.96*sd/sqrt(n).  The two-sample form would be far too lax.

SPIN-UP.  At the 16-year preliminary four years were discarded, not the campaign's
usual two: RSBLB changes the
surface energy balance and deep soil temperature has multi-year memory.  The trend column
is printed so a still-drifting difference is visible rather than averaged away.

VERIFIED BEFORE SCORING (the W1 lesson).  The staged libarpifs.SP.so in each run dir is
the 2026-08-21 15:33 rebuild and carries RSBLB; NODE.001_01 echoes
  SB0 iostat=-1 RSBLB=5.0   SB1 iostat=0 RSBLB=3.0   SB2 iostat=0 RSBLB=2.0
so the group is genuinely absent in the control and genuinely read in the arms.
NODE.001_01 is ~1.9 GB and contains NUL bytes -- grep it with -a or it silently reports
nothing and looks like a failed gate.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import sys
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')
from runs import RT, LSMF

CTL = ('SB0', 'amip_SB0_ctl')
ARMS = [('SB1 b=3', 'amip_SB1_sblb3'), ('SB2 b=2', 'amip_SB2_sblb2'),
        ('SB3 lmin120', 'amip_SB3_lmin120')]
# Campaign-standard window once the arms are complete; the 16-year preliminary used
# 1874-1885.  Override on the command line as: sb_series_eval.py Y0 Y1
Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1872, 1915)
DJF, JJA = [11, 0, 1], [5, 6, 7]

with xr.open_dataset(LSMF, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values)
    m = m[0] if m.ndim == 3 else m
    lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
land = m > 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], m.shape).copy()
arctic = land & np.broadcast_to(((lat >= 60) & (lat <= 90))[:, None], m.shape)
mid = land & np.broadcast_to(((lat >= 30) & (lat < 60))[:, None], m.shape)
sib = land & np.broadcast_to(((lat >= 50) & (lat <= 70))[:, None], m.shape) \
           & np.broadcast_to(((lon >= 60) & (lon <= 140))[None, :], m.shape)
glob = np.ones_like(m, dtype=bool)
trop = np.broadcast_to(((lat >= -30) & (lat < 30))[:, None], m.shape)
so = np.broadcast_to(((lat >= -65) & (lat < -45))[:, None], m.shape)


def am(f, s):
    k = s & np.isfinite(f)
    return float(np.average(f[k], weights=W[k]))


def load(run, var, plev=None, tag=''):
    out = []
    for y in range(Y0, Y1 + 1):
        p = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m{tag}_{y}-{y}.nc'
        if not os.path.exists(p):
            return None
        with xr.open_dataset(p, decode_times=False) as d:
            k = 't' if plev is not None else \
                [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = d[k]
            if plev is not None:
                pl = np.squeeze(d['pressure_levels'].values)
                dim = [x for x in a.dims if 'press' in x][0]
                a = a.isel({dim: int(np.argmin(np.abs(pl - plev)))})
            out.append(np.squeeze(a.values))
    return np.stack(out)


def stat(dy):
    """mean, 95% paired threshold, and the per-decade trend of the difference."""
    n = np.isfinite(dy).sum()
    mu = float(np.nanmean(dy))
    thr = 1.96 * float(np.nanstd(dy, ddof=1)) / np.sqrt(n)
    x = np.arange(len(dy))
    slope = float(np.polyfit(x[np.isfinite(dy)], dy[np.isfinite(dy)], 1)[0]) * 10
    return mu, thr, slope


C = {v: load(CTL[1], *v) for v in [('2t',), ('skt',), ('pl_t', 92500, '_pl'),
                                   ('pl_t', 85000, '_pl'), ('tsr',), ('ttr',),
                                   ('tsrc',)]}
if any(v is None for v in C.values()):
    sys.exit('control incomplete')

print(f'SB series vs {CTL[0]}, AMIP paired, DJF {Y0}-{Y1} ({Y1-Y0+1} yr, 4 discarded).')
print('* = beyond the paired 95 % threshold.  trend = K/decade in the DIFFERENCE.\n')

for label, run in ARMS:
    A = {v: load(run, *v) for v in C}
    if any(v is None for v in A.values()):
        print(f'{label}: incomplete\n')
        continue
    print(f'=== {label} ===')
    print(f'{"metric":>28} | {"delta":>8} {"thr":>6} {"trend":>7}')
    rows = []
    for nm, sel in [('60-90N land', arctic), ('30-60N land', mid), ('all land', land)]:
        # the primary metric: change in the low-level inversion
        dy = np.array([am((A[('pl_t',92500,'_pl')][i][DJF].mean(0) - A[('2t',)][i][DJF].mean(0))
                        - (C[('pl_t',92500,'_pl')][i][DJF].mean(0) - C[('2t',)][i][DJF].mean(0)), sel)
                       for i in range(len(A[('2t',)]))])
        rows.append((f'DJF inversion {nm}', *stat(dy)))
    for nm, sel in [('60-90N land', arctic), ('30-60N land', mid), ('all land', land)]:
        dy = np.array([am(A[('2t',)][i][DJF].mean(0) - C[('2t',)][i][DJF].mean(0), sel)
                       for i in range(len(A[('2t',)]))])
        rows.append((f'DJF T2m {nm}', *stat(dy)))
    # the layer above: a screen warming paid for by cooling aloft is the signature we want
    dy = np.array([am(A[('pl_t',92500,'_pl')][i][DJF].mean(0) - C[('pl_t',92500,'_pl')][i][DJF].mean(0), arctic)
                   for i in range(len(A[('2t',)]))])
    rows.append(('DJF T925 60-90N land', *stat(dy)))
    dy = np.array([am(A[('pl_t',85000,'_pl')][i][DJF].mean(0) - C[('pl_t',85000,'_pl')][i][DJF].mean(0), arctic)
                   for i in range(len(A[('2t',)]))])
    rows.append(('DJF T850 60-90N land', *stat(dy)))
    # skin: should barely move, it was already right
    dy = np.array([am(A[('skt',)][i][DJF].mean(0) - C[('skt',)][i][DJF].mean(0), arctic)
                   for i in range(len(A[('2t',)]))])
    rows.append(('DJF skt 60-90N land', *stat(dy)))
    # guardrails
    dy = np.array([am(A[('2t',)][i][JJA].mean(0) - C[('2t',)][i][JJA].mean(0), sib)
                   for i in range(len(A[('2t',)]))])
    rows.append(('GUARD Siberia JJA T2m', *stat(dy)))
    dy = np.array([am((A[('tsr',)][i] + A[('ttr',)][i]).mean(0) / 3600.0
                    - (C[('tsr',)][i] + C[('ttr',)][i]).mean(0) / 3600.0, glob)
                   for i in range(len(A[('2t',)]))])
    rows.append(('GUARD global net TOA', *stat(dy)))
    # The tropical shortwave cost is why SB3 exists: SB2 pays +0.515 of SW CRE there,
    # concentrated in the subtropical Sc decks, which lands on nino34.
    for rn, sel in [('tropics', trop), ('SO 45-65S', so)]:
        dy = np.array([am(((A[('tsr',)][i] - A[('tsrc',)][i])
                         - (C[('tsr',)][i] - C[('tsrc',)][i])).mean(0) / 3600.0, sel)
                       for i in range(len(A[('2t',)]))])
        rows.append((f'COST {rn} SW CRE', *stat(dy)))
    dy = np.array([am((A[('tsr',)][i] + A[('ttr',)][i]).mean(0) / 3600.0
                    - (C[('tsr',)][i] + C[('ttr',)][i]).mean(0) / 3600.0, trop)
                   for i in range(len(A[('2t',)]))])
    rows.append(('COST tropics net TOA', *stat(dy)))
    for nm, mu, thr, sl in rows:
        star = '*' if abs(mu) > thr else ' '
        print(f'{nm:>28} | {mu:+8.3f}{star} {thr:6.3f} {sl:+7.3f}')
    print()
