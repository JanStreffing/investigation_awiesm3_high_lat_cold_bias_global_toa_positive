"""Atmospheric freshwater over the Southern Ocean and Antarctica from OpenIFS output.

Annual means over Y0-Y1 of lsp+cp (precipitation), e (evaporation, negative = upward) and
ro (runoff) on the remapped grid, accumulated m per hourly step (/3600 -> m/s), cos-lat
weighted.  Integrated over the ocean 60-78S (P, E, P-E) and over Antarctic land south of
60S (P-E and ro, the mass that the coupling returns to the ocean as calving/runoff), in
Gt/yr and mm/yr.  The land-sea mask comes from PI200's lsm file (same TCO95 grid).

Usage:  ROOT=... ARM=PI Y0=2000 Y1=2014 python3 scripts/analysis/so_atm_freshwater_nodes.py
"""
import os, glob
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ['ROOT']; ARM = os.environ['ARM']; Y0, Y1 = int(os.environ['Y0']), int(os.environ['Y1'])
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'


def ann(v):
    acc = 0; n = 0
    for y in range(Y0, Y1 + 1):
        p = glob.glob(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc') + glob.glob(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc')
        if not p: continue
        with xr.open_dataset(p[0], decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8'); lat = np.squeeze(d['lat'].values)
        acc = acc + a.mean(0) / float(os.environ.get('ACC', 3600.0)); n += 1     # ACC=1 for runs whose output is already a rate (awicm3 v3.3)
    return (acc / n if n else None), lat


P, lat = ann('lsp'); C, _ = ann('cp'); E, _ = ann('e'); RO, _ = ann('ro')
P = P + C
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = lsm[0] if lsm.ndim == 3 else lsm
dlat = np.abs(np.gradient(lat)); A = (np.cos(np.deg2rad(lat)) * dlat * 111.195e3 * (360.0 / P.shape[1]) * 111.195e3)[:, None] * np.ones(P.shape[1])[None, :]
GT = 1000.0 * 365.25 * 86400 / 1e12
ocn = (lsm <= 0.5) & (lat[:, None] >= -78) & (lat[:, None] <= -60) & np.ones_like(P, bool)
lnd = (lsm > 0.5) & (lat[:, None] < -60) & np.ones_like(P, bool)
mm = lambda f, k: np.sum(f[k] * A[k]) / A[k].sum() * 365.25 * 86400 * 1000
print(f'{ARM} {Y0}-{Y1} (OpenIFS): ocean 60-78S {A[ocn].sum() / 1e12:.2f} M km2, Antarctic land {A[lnd].sum() / 1e12:.2f} M km2')
print(f'  ocean 60-78S: P {np.sum(P[ocn] * A[ocn]) * GT:6.0f} Gt/yr ({mm(P, ocn):.0f} mm/yr), E {np.sum(E[ocn] * A[ocn]) * GT:6.0f} ({mm(E, ocn):.0f}), P+E {np.sum((P + E)[ocn] * A[ocn]) * GT:6.0f} ({mm(P + E, ocn):.0f})')
print(f'  Antarctic land: P+E {np.sum((P + E)[lnd] * A[lnd]) * GT:6.0f} Gt/yr, runoff ro {np.sum(RO[lnd] * A[lnd]) * GT:6.0f} Gt/yr' if RO is not None else '  Antarctic land: no ro output')
