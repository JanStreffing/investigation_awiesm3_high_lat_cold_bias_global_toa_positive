"""How many AMIP years does the SB screen actually need?

WHY.  SB0/SB1/SB2 are 46-year arms in 8-year legs.  46 years is the campaign default, not
a derived number, and the decision the screen has to support -- "does RSBLB move the
winter inversion, and does it break anything" -- may be answerable much sooner, or may
need the full length for the guardrails rather than for the signal.

METHOD.  The noise is measured, not assumed: take an existing AMIP pair of the same
construction (J1 - I1, 44 yr, one namelist number apart) and compute the standard
deviation of the ANNUAL paired difference for each metric.  The paired 95 % detection
threshold at n years is then 1.96 * sd / sqrt(n).  This is the right form for AMIP arms
under identical prescribed SST; the two-sample 1.96*sd*sqrt(2/n) would be far too lax.

The target sizes come from amip_vs_coupled_inversion.py: the AMIP low-level inversion
bias is +1.54 K at 60-90N, so a lever removing half of it must resolve 0.77 K and one
removing a quarter must resolve 0.39 K.

RESULT (2026-08-21).  SCORE THE INVERSION, NOT THE SCREEN TEMPERATURE.  At 60-90N the
paired annual sd is 1.377 K on T2m but only 0.611 K on T925-T2m, because differencing two
levels in the same column cancels the synoptic variability that dominates Arctic winter
T2m.  That is a factor 2.25 in noise and a factor ~5 in years:

    detect                  on T2m     on the inversion
    half the bias (0.77 K)   12.3 yr        2.4 yr
    a quarter    (0.39 K)    49.1 yr        9.7 yr

So the go/no-go on RSBLB is available at leg 2 (~16 yr, threshold 0.32 K on the
inversion).  What needs the full 46 years is the SIBERIA JJA GUARDRAIL -- 0.40 K at n=14
against 0.23 K at n=44 -- not the signal.  Global net TOA is inside its +-0.3 guardrail
by leg 1.

And do not serialise the coupled pair behind AMIP finishing: the AMIP-to-coupled transfer
ratio in this campaign has ranged 0.19x to 1.01x by lever, so extra AMIP precision buys
sign and guardrails, never a coupled magnitude.

At a short read, discard 3-4 spin-up years rather than the campaign's usual 2: a mixing
change moves the surface energy balance and deep soil temperature has multi-year memory.
Check the inversion difference is not still trending before calling it.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')
from runs import RT, LSMF, Y0, Y1

CTL, ARM = 'amip_I1_scf', 'amip_J1_lamsk15'
SEAS = {'DJF': [11, 0, 1], 'JJA': [5, 6, 7]}

with xr.open_dataset(LSMF, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values)
    m = m[0] if m.ndim == 3 else m
    lat = np.squeeze(d['lat'].values)
    lon = np.squeeze(d['lon'].values)
land = m > 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], m.shape).copy()
# Siberia as the campaign defines it
sib = land & np.broadcast_to(((lat >= 50) & (lat <= 70))[:, None], m.shape) \
           & np.broadcast_to(((lon >= 60) & (lon <= 140))[None, :], m.shape)
arctic = land & np.broadcast_to(((lat >= 60) & (lat <= 90))[:, None], m.shape)


def am(f, s):
    k = s & np.isfinite(f)
    return float(np.average(f[k], weights=W[k]))


def series(run, var):
    out = []
    for y in range(Y0, Y1 + 1):
        p = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        with xr.open_dataset(p, decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            out.append(np.squeeze(d[k].values))
    return np.stack(out)


t_c, t_a = series(CTL, '2t'), series(ARM, '2t')
# TOA net: accumulated fluxes, /3600 for monthly means (verified: global ASR 240.4)
tsr_c, tsr_a = series(CTL, 'tsr'), series(ARM, 'tsr')
ttr_c, ttr_a = series(CTL, 'ttr'), series(ARM, 'ttr')

METRICS = {}
for nm, sel, mons in [('60-90N land DJF T2m', arctic, SEAS['DJF']),
                      ('Siberia JJA T2m', sib, SEAS['JJA']),
                      ('all-land DJF T2m', land, SEAS['DJF'])]:
    METRICS[nm] = np.array([am(t_a[i][mons].mean(0) - t_c[i][mons].mean(0), sel)
                            for i in range(t_a.shape[0])])
glob = np.ones_like(m, dtype=bool)
METRICS['global net TOA'] = np.array(
    [am((tsr_a[i] + ttr_a[i]).mean(0) / 3600.0 - (tsr_c[i] + ttr_c[i]).mean(0) / 3600.0, glob)
     for i in range(t_a.shape[0])])

print(f'Paired annual-difference noise, measured on {ARM} - {CTL}, {Y0}-{Y1}\n')
print(f'{"metric":>22} | {"sd":>6} | ' + ' | '.join(f'n={n:<3}' for n in (8, 14, 22, 30, 44)))
print('-' * 78)
for nm, d in METRICS.items():
    sd = float(np.nanstd(d, ddof=1))
    thr = [1.96 * sd / np.sqrt(n) for n in (8, 14, 22, 30, 44)]
    print(f'{nm:>22} | {sd:6.3f} | ' + ' | '.join(f'{t:5.3f}' for t in thr))

print('\nYears needed to resolve a given response at 60-90N land DJF T2m:')
sd = float(np.nanstd(METRICS['60-90N land DJF T2m'], ddof=1))
for frac, label in [(1.00, 'all of the +1.54 K AMIP inversion bias'),
                    (0.50, 'half of it (0.77 K)'),
                    (0.25, 'a quarter of it (0.39 K)'),
                    (0.10, 'a tenth of it (0.15 K)')]:
    target = 1.54 * frac
    n = (1.96 * sd / target) ** 2
    print(f'  {target:5.2f} K  ->  n = {n:6.1f} yr   ({label})')
