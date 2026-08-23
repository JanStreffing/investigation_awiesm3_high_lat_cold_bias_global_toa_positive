"""11Q/11R: does the AMIP stable-BL trade-off survive coupling?

THE QUESTION.  SB2 (RSBLB 2.0) in AMIP gave a clean but bounded result: Arctic DJF
inversion -0.713 K, mid-latitude land T2m +0.492 K, Siberian JJA untouched -- but the
Arctic SCREEN did not move significantly (+0.285, thr 0.345), and it billed +0.175 W/m2
of global net TOA, essentially all of it tropical shortwave from low-cloud loss
(tropical SW CRE +0.515).  SB3 then showed the exchange rate is a property of the scheme
rather than the knob: tropical SW CRE per K of Arctic inversion is 0.88 (SB1), 0.87 (SB2),
3.43 (SB3, the mixing-length route).  Constant across an RSBLB dose change, so it cannot
be dodged by tuning gentler.

Two things could only be settled coupled:
  1. Does the Arctic screen follow once SEA ICE is free to respond?  AMIP prescribes ice,
     which suppresses the one feedback that could amplify an Arctic surface warming, and
     45 % of the coupled Arctic inversion bias has no AMIP counterpart at all
     (amip_vs_coupled_inversion.py: +1.54 K AMIP against +2.81 K coupled).
  2. What does the tropical shortwave cost DO once the ocean can respond?  In AMIP it is
     an inert imbalance; coupled it warms the tropical ocean and moves nino34, which is
     already where cmpitool scored 11G worse than 11E.

PAIRS.  11Q - 11N (1850 forcing) and 11R - 11P (1990).  Each arm is its control plus one
namelist number, branched from the same 1350 state, on a binary whose only difference is
the &NAMVDFS patch -- bit-identical with the group absent, and the coupled source has no
commits since the 08-12 build the controls used.  So these are true one-parameter pairs
and the statistics are paired.

SEA ICE IS SCORED HERE because it is the guardrail AMIP could not test even in principle.

CAVEAT ON 11Q.  Its leg 1380-89 was lost once to a transient Lustre read failure
(LPJ-GUESS "failed to find element to deserialize" on a state file that reads cleanly
afterwards; DKRZ Infiniband/Lustre incident open since 2026-08-17) and was re-run.  The
failed attempt is preserved at run_13800101-13891231_dead_lustre_20260822.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import sys
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
PAIRS = [('11Q - 11N  (1850)', f'{R}/11Q', f'{R}/11N'),
         ('11R - 11P  (1990)', f'{R}/11R', f'{R}/11P')]
Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1360, 1389)
DJF, JJA, SON = [11, 0, 1], [5, 6, 7], [8, 9, 10]

f = sorted(glob.glob(f'{R}/11P/outdata/oifs/atm_remapped_1m_lsm_*.nc'))[0]
with xr.open_dataset(f, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values)
    m = m[0] if m.ndim == 3 else m
    lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
land, ocean = m > 0.5, m <= 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], m.shape).copy()
bnd = lambda lo, hi: np.broadcast_to(((lat >= lo) & (lat < hi))[:, None], m.shape)
arctic, mid = land & bnd(60, 90), land & bnd(30, 60)
sib = land & bnd(50, 70) & np.broadcast_to(((lon >= 60) & (lon <= 140))[None, :], m.shape)
trop, so = bnd(-30, 30), bnd(-65, -45)
nh_ice, sh_ice = ocean & bnd(45, 90), ocean & bnd(-90, -45)
glob_ = np.ones_like(m, dtype=bool)
# grid-cell area for ice AREA rather than mean concentration
RE = 6.371e6
dlat = np.deg2rad(abs(lat[1] - lat[0])); dlon = 2 * np.pi / m.shape[1]
AREA = np.broadcast_to((RE ** 2 * np.cos(np.deg2rad(lat)) * dlat * dlon)[:, None], m.shape)


def am(f_, s):
    k = s & np.isfinite(f_)
    return float(np.average(f_[k], weights=W[k])) if k.any() else np.nan


def load(path, var, plev=None):
    out = []
    for y in range(Y0, Y1 + 1):
        p = f'{path}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
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
    dy = np.asarray(dy, dtype=float)
    n = np.isfinite(dy).sum()
    mu = float(np.nanmean(dy))
    thr = 1.96 * float(np.nanstd(dy, ddof=1)) / np.sqrt(n)
    x = np.arange(len(dy))
    sl = float(np.polyfit(x, dy, 1)[0]) * 10
    return mu, thr, sl


VARS = [('2t', None), ('skt', None), ('pl_t', 92500), ('pl_t', 85000),
        ('tsr', None), ('tsrc', None), ('ttr', None), ('ci', None), ('sp', None)]

print(f'Coupled SB pairs, {Y0}-{Y1} ({Y1-Y0+1} yr), paired.  '
      f'* = beyond the paired 95 % threshold.\n')

for label, ap, cp in PAIRS:
    A = {v: load(ap, *v) for v in VARS}
    C = {v: load(cp, *v) for v in VARS}
    if any(v is None for v in A.values()) or any(v is None for v in C.values()):
        miss = [v for v in VARS if A[v] is None or C[v] is None]
        print(f'{label}: incomplete ({miss})\n')
        continue
    n = A[('2t', None)].shape[0]
    sp = np.mean([C[('sp', None)][i][DJF].mean(0) for i in range(n)], axis=0)
    rows = []

    def add(nm, dy):
        rows.append((nm, *stat(dy)))

    # --- the target: does the inversion weaken, and does the SCREEN follow? ---
    for rn, sel in [('60-90N', arctic), ('30-60N', mid), ('all land', land)]:
        s9 = sel & (sp >= 92500)
        add(f'DJF inversion {rn}', [am((A[('pl_t',92500)][i][DJF].mean(0) - A[('2t',None)][i][DJF].mean(0))
                                     - (C[('pl_t',92500)][i][DJF].mean(0) - C[('2t',None)][i][DJF].mean(0)), s9)
                                    for i in range(n)])
    for rn, sel in [('60-90N', arctic), ('30-60N', mid), ('all land', land)]:
        add(f'DJF T2m {rn}', [am(A[('2t',None)][i][DJF].mean(0) - C[('2t',None)][i][DJF].mean(0), sel)
                              for i in range(n)])
    add('DJF T925 60-90N', [am(A[('pl_t',92500)][i][DJF].mean(0) - C[('pl_t',92500)][i][DJF].mean(0),
                               arctic & (sp >= 92500)) for i in range(n)])
    add('DJF T850 60-90N', [am(A[('pl_t',85000)][i][DJF].mean(0) - C[('pl_t',85000)][i][DJF].mean(0),
                               arctic & (sp >= 85000)) for i in range(n)])
    add('ANN T2m all land', [am(A[('2t',None)][i].mean(0) - C[('2t',None)][i].mean(0), land) for i in range(n)])
    add('ANN T2m global', [am(A[('2t',None)][i].mean(0) - C[('2t',None)][i].mean(0), glob_) for i in range(n)])
    # --- guardrails ---
    add('GUARD Siberia JJA', [am(A[('2t',None)][i][JJA].mean(0) - C[('2t',None)][i][JJA].mean(0), sib)
                              for i in range(n)])
    add('GUARD global net TOA', [am((A[('tsr',None)][i]+A[('ttr',None)][i]).mean(0)/3600.
                                  - (C[('tsr',None)][i]+C[('ttr',None)][i]).mean(0)/3600., glob_)
                                 for i in range(n)])
    for rn, sel in [('tropics', trop), ('SO 45-65S', so)]:
        add(f'COST {rn} SW CRE', [am(((A[('tsr',None)][i]-A[('tsrc',None)][i])
                                    - (C[('tsr',None)][i]-C[('tsrc',None)][i])).mean(0)/3600., sel)
                                  for i in range(n)])
    add('COST tropics net TOA', [am((A[('tsr',None)][i]+A[('ttr',None)][i]).mean(0)/3600.
                                  - (C[('tsr',None)][i]+C[('ttr',None)][i]).mean(0)/3600., trop)
                                 for i in range(n)])
    add('COST tropics SST-proxy', [am(A[('skt',None)][i].mean(0) - C[('skt',None)][i].mean(0), trop & ocean)
                                   for i in range(n)])
    # --- sea ice: the guardrail AMIP cannot test ---
    for rn, sel, mons in [('NH ice Mar', nh_ice, [2]), ('NH ice Sep', nh_ice, [8]),
                          ('SH ice Sep', sh_ice, [8]), ('SH ice Mar', sh_ice, [2])]:
        add(f'ICE {rn} [1e6 km2]',
            [((A[('ci',None)][i][mons].mean(0) - C[('ci',None)][i][mons].mean(0)) * AREA)[sel].sum() / 1e12
             for i in range(n)])
    print(f'=== {label} ===')
    print(f'{"metric":>26} | {"delta":>9} {"thr":>7} {"trend/dec":>10}')
    for nm, mu, thr, sl in rows:
        print(f'{nm:>26} | {mu:+9.3f}{"*" if abs(mu) > thr else " "} {thr:7.3f} {sl:+10.3f}')
    print()
