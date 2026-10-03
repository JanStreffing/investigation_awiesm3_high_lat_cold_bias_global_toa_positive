"""Did the GM change reach the run?  Bolus velocities of the branches against the control.

Annual-mean magnitude of the GM bolus velocity (bolus_u, bolus_v, bolus_w on FESOM's regular
output) in 100-500 m, area-weighted by latitude band, branch over control.  With the
Rossby cutoff replacing the resolution scaling the coefficient rises about 2.7x at 60-78S
and 2.5-3.5x north of 45N, and halves in the tropics; the bolus velocity scales with it.

Usage:  CONTROL=PI200 RUNS=PI200_gmR2500,PI200_gmR1500 YEAR=1580 python3 scripts/analysis/gm_bolus_check.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
CONTROL = os.environ.get('CONTROL', 'PI200'); RUNS = os.environ.get('RUNS', 'PI200_gmR2500,PI200_gmR1500').split(',')
Y = int(os.environ.get('YEAR', 1580))
BANDS = ((-78, -60), (-60, -45), (-45, -30), (-30, 30), (30, 45), (45, 60), (60, 80))


def load(arm):
    out = {}
    for v in ('bolus_u', 'bolus_v', 'bolus_w'):
        with xr.open_dataset(f'{R}/{arm}/outdata/fesom/{v}.fesom.gr.{Y}.nc', decode_times=False) as d:
            a = d[v].mean('time').values.astype('f8'); lat = d['lat'].values
            zz = d['nz'].values if 'nz' in d else d['nz1'].values          # bolus_w sits on the 48 interfaces, u/v on the 47 mid-levels
        a[np.abs(a) > 100] = np.nan; out[v] = a[..., (zz >= 100) & (zz < 500)]
    return out, lat


C, lat = load(CONTROL)
W = np.cos(np.deg2rad(lat))[:, None] * np.ones(C['bolus_u'].shape[1])[None, :]
mag = lambda F: np.nanmean(np.sqrt(F['bolus_u'] ** 2 + F['bolus_v'] ** 2), -1)
wmag = lambda F: np.nanmean(np.abs(F['bolus_w']), -1)
cm, cw = mag(C), wmag(C)
print(f'{Y}: |bolus horizontal| and |bolus_w| in 100-500 m, branch / control ({CONTROL}), by band')
print(f'  {"band":<10}{"ctl |u| mm/s":>13}' + ''.join(f'{a.replace(CONTROL + "_", ""):>12}' for a in RUNS) + f'{"ctl |w| um/s":>14}' + ''.join(f'{a.replace(CONTROL + "_", ""):>12}' for a in RUNS))
for lo, hi in BANDS:
    k = (lat[:, None] >= lo) & (lat[:, None] < hi) & np.isfinite(cm) & (cm > 0)
    row = f'  {lo:>4}..{hi:<4}{np.average(cm[k], weights=W[k]) * 1e3:13.3f}'
    ratios_w = ''
    for a in RUNS:
        F, _ = load(a); bm, bw = mag(F), wmag(F)
        row += f'{np.average(bm[k], weights=W[k]) / np.average(cm[k], weights=W[k]):12.2f}'
        ratios_w += f'{np.average(bw[k], weights=W[k]) / np.average(cw[k], weights=W[k]):12.2f}'
    print(row + f'{np.average(cw[k], weights=W[k]) * 1e6:14.3f}' + ratios_w)
