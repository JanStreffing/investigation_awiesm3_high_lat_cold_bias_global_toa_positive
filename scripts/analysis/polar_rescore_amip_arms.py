"""Task 1: re-score every AMIP cloud arm on the POLAR OCEAN band it was never read on.

WHY.  The corrected target is +13.36 W/m2 of DJF SW CRE over 90-60S OCEAN, and roughly
HALF of it is cloud AMOUNT rather than opacity (cover 5.8 pp low, each cloud 6.7 % dim).
The campaign owns several amount levers -- B2 (RCLDIFF_CONVI), B3 (RCLDIFF), A1b
(RCL_OVERLAPLIQICE), B4 (ENTSHALP) -- and every one was judged on Siberian JJA, global
surface flux, or the 45-65S SO metric.  None was ever read on the polar ocean.

That is exactly how A1c came to be recorded as "failed -- cloud-depth selectivity idea
dead" when re-scoring showed it significant in the target band (-0.815 DJF), null in the
Arctic (+0.017 against a 0.593 threshold) and free on TOA (+0.006).  The verdicts were
rendered against the wrong band, not wrongly computed.

WHAT DECIDES USABILITY, in order:
  1. does it brighten 90-60S ocean in DJF (the target)
  2. does it leave 60-90N ocean JJA alone (the Arctic is already -18 overcorrected, and
     this killed SC3)
  3. does it avoid driving 60-45S ocean further negative (already -3.70)
  4. tundra and Siberia untouched, global net TOA within guardrail

All bands are OCEAN-MASKED: averaging whole latitude bands folds the Antarctic ice sheet
into 90-60S and halves the apparent signal.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')
from runs import RT, LSMF, OBS, Y0, Y1

CTL = 'amip_pi_base'
ARMS = [('A1a ovl0.10', 'amip_A1_overlap01'), ('A1b ovl0.35', 'amip_A1_overlap035'),
        ('A1c depth1500', 'amip_A1c_depliqdepth1500'), ('A2 KKland150', 'amip_A2_kknumland150'),
        ('B2 convi25', 'amip_B2_clddiffconvi25'), ('B3 clddiff', 'amip_B3_clddiff15e6'),
        ('B4 entshalp3', 'amip_B4_entshalp3'), ('B5 capdcycl0', 'amip_B5_capdcycl0'),
        ('B6 lcritsnow', 'amip_B6_lcritsnow1e5'), ('B7 rvice0.22', 'amip_B7_rvice022')]

with xr.open_dataset(LSMF, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values)
    m = m[0] if m.ndim == 3 else m
    lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
ocean, land = m <= 0.5, m > 0.5
RE = 6.371e6
A = np.broadcast_to((RE**2 * np.cos(np.deg2rad(lat)) * np.deg2rad(abs(lat[1] - lat[0]))
                     * 2 * np.pi / m.shape[1])[:, None], m.shape)
sib = land & np.broadcast_to(((lat >= 50) & (lat <= 70))[:, None], m.shape) \
           & np.broadcast_to(((lon >= 60) & (lon <= 140))[None, :], m.shape)
DJF, JJA, ANN = [11, 0, 1], [5, 6, 7], list(range(12))


def yrs(run, var):
    out = []
    for y in range(Y0, Y1 + 1):
        p = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(p):
            return None
        with xr.open_dataset(p, decode_times=False) as d:
            out.append(np.squeeze(d[[k for k in d.data_vars
                                     if 'bnds' not in k and 'bounds' not in k][0]].values))
    return np.stack(out)


def am(f, s):
    k = s & np.isfinite(f)
    return float(np.average(f[k], weights=A[k]))


def stat(dy):
    dy = np.asarray(dy)
    return float(np.nanmean(dy)), 1.96 * float(np.nanstd(dy, ddof=1)) / np.sqrt(len(dy))


V = ('tsr', 'tsrc', 'ttr', '2t')
C = {v: yrs(CTL, v) for v in V}
BANDS = [('polar 90-60S DJF', ocean & np.broadcast_to((lat < -60)[:, None], m.shape), DJF),
         ('ARCTIC 60-90N JJA', ocean & np.broadcast_to((lat > 60)[:, None], m.shape), JJA),
         ('subAnt 60-45S DJF', ocean & np.broadcast_to(((lat >= -60) & (lat < -45))[:, None], m.shape), DJF),
         ('tropics ocean ANN', ocean & np.broadcast_to(((lat >= -30) & (lat < 30))[:, None], m.shape), ANN)]

print(__doc__)
print('=' * 104)
print(f'AMIP arms minus control, OCEAN-masked, paired, {Y0}-{Y1}.  '
      f'SW CRE in W/m2; * beyond 95 % threshold.\n')
hdr = f'{"arm":>15} |'
for b, _, _ in BANDS:
    hdr += f' {b:>18}'
hdr += f' | {"SibJJA":>8} {"netTOA":>8}'
print(hdr)
print('-' * 104)
for lab, run in ARMS:
    Aa = {v: yrs(run, v) for v in V}
    if any(x is None for x in Aa.values()):
        print(f'{lab:>15} | incomplete'); continue
    n = min(C['tsr'].shape[0], Aa['tsr'].shape[0])
    line = f'{lab:>15} |'
    for _, sel, mo in BANDS:
        mu, th = stat([am(((Aa['tsr'][i][mo].mean(0) - Aa['tsrc'][i][mo].mean(0))
                          - (C['tsr'][i][mo].mean(0) - C['tsrc'][i][mo].mean(0))) / 3600., sel)
                       for i in range(n)])
        line += f' {mu:+13.2f}{"*" if abs(mu) > th else " "}({th:4.2f})'
    mu, th = stat([am(Aa['2t'][i][JJA].mean(0) - C['2t'][i][JJA].mean(0), sib) for i in range(n)])
    line += f' | {mu:+7.2f}{"*" if abs(mu) > th else " "}'
    g = np.ones_like(m, dtype=bool)
    mu, th = stat([am(((Aa['tsr'][i] + Aa['ttr'][i]).mean(0)
                      - (C['tsr'][i] + C['ttr'][i]).mean(0)) / 3600., g) for i in range(n)])
    line += f' {mu:+7.2f}{"*" if abs(mu) > th else " "}'
    print(line)
print('\nTARGET: polar 90-60S DJF must go NEGATIVE (model reflects too little, +13.36).')
print('VETO:   ARCTIC 60-90N JJA must stay NULL (already -18.22 overcorrected).')
print('VETO:   subAnt 60-45S DJF must not go further negative (already -3.70).')
