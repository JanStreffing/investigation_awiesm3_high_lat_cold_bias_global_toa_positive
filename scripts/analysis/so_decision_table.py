"""Which configuration puts the Southern Ocean SURFACE closest to zero for the least
collateral damage?

THE QUESTION, stated so it can be scored.  Three constraints, in the order they bind:

  1. SO 45-65S OCEAN SURFACE net SW, bias vs CERES -> as close to ZERO as possible.
     This is the field the campaign's "SW RMSE" is built on and the one the coupled
     ocean actually feels.  It is NOT the TOA cloud forcing that most of this campaign
     has been quoting; the two disagree, and both are reported here so the disagreement
     stays visible.
  2. TROPICS -- do not spend them.  The tropics sit BELOW CERES (control net TOA 42.61
     against 45.11, -0.67 period-clean), so they need MORE absorbed energy, not less.
     The tolerance is therefore ONE-SIDED: a NEGATIVE dTrop is a cost, a positive one is
     an improvement.  An earlier version of this table flagged |dTrop|>0.5 regardless of
     sign and so marked LX1 as failing for moving the tropics the RIGHT way.
  3. SIBERIA -- do not give back the land campaign's gains.  JJA is what the forest
     needs (the GDD gate); the 44-yr thresholds are +-0.246 JJA and +-0.604 DJF.

Global net TOA is carried alongside because a configuration that nails the Southern
Ocean by dimming the planet is not a Southern Ocean fix -- LX3 is exactly that case.

ALL 44-YEAR ARMS with output are scored, not just the ones I expected to win, because
the point of the exercise is to let an unexpected arm surface.  S4 in particular was
completed on 2026-08-09 and never scored until round 27.

WINDOW 1872-1915, CERES EBAF climatology, ocean points only for the SO surface metric.
Deltas for tropics and Siberia are against amip_pi_base; the SO and surface numbers are
ABSOLUTE biases vs CERES, because "closest to zero" is a statement about the bias, not
about the change.
"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runs import RUNS, RT, Y0, Y1, LSMF, OBS

ACC = 3600.0
SO = (-65.0, -45.0)
SIB = (55.0, 75.0, 60.0, 180.0)
TOL_TROP, THR_JJA, THR_DJF = 0.50, 0.246, 0.604

print(__doc__)
print('=' * 118)

with xr.open_dataset(OBS) as c:
    ce_sfc = c['sfc_net_sw_all_clim'].values.mean(axis=0)
    ce_cre = c['toa_cre_sw_clim'].values.mean(axis=0)
    clat, clon = c['lat'].values, c['lon'].values
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values)
    if lsm.ndim == 3:
        lsm = lsm[0]
    llat, llon = d['lat'].values, d['lon'].values


def load(run, var):
    acc, lat, lon = [], None, None
    for y in range(Y0, Y1 + 1):
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            return None, None, None
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values
            lat, lon = d['lat'].values, d['lon'].values
        if a.shape[0] != 12:
            return None, None, None
        acc.append(a / (ACC if var in ('ssr', 'tsr', 'tsrc', 'ttr', 'ttrc') else 1.0))
    return np.mean(acc, axis=0), lat, lon


def score(run):
    ssr, lat, lon = load(run, 'ssr')
    if ssr is None:
        return None
    tsr, _, _ = load(run, 'tsr'); tsrc, _, _ = load(run, 'tsrc')
    ttr, _, _ = load(run, 'ttr'); t2, _, _ = load(run, '2t')
    if any(x is None for x in (tsr, tsrc, ttr, t2)):
        return None
    ii = np.abs(clat[None, :] - lat[:, None]).argmin(axis=1)
    jj = np.abs(clon[None, :] - lon[:, None]).argmin(axis=1)
    li = np.abs(llat[None, :] - lat[:, None]).argmin(axis=1)
    lj = np.abs(llon[None, :] - lon[:, None]).argmin(axis=1)
    land = lsm[np.ix_(li, lj)] > 0.5

    so = (lat >= SO[0]) & (lat < SO[1])
    W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], ssr.shape[1:])

    b_sfc = ssr.mean(axis=0) - ce_sfc[np.ix_(ii, jj)]
    m = so[:, None] & (~land)
    so_sfc = float(np.average(b_sfc[m], weights=W[m]))
    so_rms = float(np.sqrt(np.average(b_sfc[m] ** 2, weights=W[m])))

    b_cre = (tsr - tsrc).mean(axis=0) - ce_cre[np.ix_(ii, jj)]
    so_cre = float(np.average(b_cre[so], weights=W[so]))

    net = (tsr + ttr).mean(axis=0)
    tr = (lat >= -30) & (lat < 30)
    trop = float(np.average(net[tr], weights=W[tr]))
    glob = float(np.average(net, weights=W))

    sy = (lat >= SIB[0]) & (lat < SIB[1]); sx = (lon >= SIB[2]) & (lon <= SIB[3])
    box = np.zeros_like(land, dtype=bool); box[np.ix_(sy, sx)] = True
    box &= land
    jja = float(np.average((t2[[5, 6, 7]].mean(axis=0) - 273.15)[box], weights=W[box]))
    djf = float(np.average((t2[[11, 0, 1]].mean(axis=0) - 273.15)[box], weights=W[box]))
    return dict(so_sfc=so_sfc, so_rms=so_rms, so_cre=so_cre, trop=trop, glob=glob,
                jja=jja, djf=djf)


res = {}
for lab, run in RUNS:
    s = score(run)
    if s:
        res[lab] = s
ctl = res.get('control')
if ctl is None:
    raise SystemExit('control missing')

rows = []
for lab, s in res.items():
    if lab == 'control':
        continue
    rows.append((lab, s['so_sfc'], s['so_rms'], s['so_cre'],
                 s['trop'] - ctl['trop'], s['jja'] - ctl['jja'], s['djf'] - ctl['djf'],
                 s['glob']))
rows.sort(key=lambda r: abs(r[1]))

print(f'CONTROL: SO surface bias {ctl["so_sfc"]:+.2f} (RMSE {ctl["so_rms"]:.2f}), '
      f'SO TOA CRE bias {ctl["so_cre"]:+.2f}, global net TOA {ctl["glob"]:+.2f}\n')
print(f'  {"run":22s} {"SO sfc":>8s} {"RMSE":>7s} {"SO CRE":>8s} | {"dTrop":>7s} '
      f'{"dJJA":>7s} {"dDJF":>7s} {"globTOA":>8s} | verdict')
print('  ' + '-' * 112)
for lab, sfc, rms, cre, dtr, djja, ddjf, gl in rows:
    bad = []
    if dtr < -TOL_TROP:                      # one-sided: only a LOSS is a cost
        bad.append(f'tropics {dtr:+.2f}')
    if djja < -THR_JJA:
        bad.append(f'JJA {djja:+.2f}')
    if ddjf < -THR_DJF:
        bad.append(f'DJF {ddjf:+.2f}')
    v = 'CLEAN' if not bad else '; '.join(bad)
    print(f'  {lab:22s} {sfc:+8.2f} {rms:7.2f} {cre:+8.2f} | {dtr:+7.2f} {djja:+7.2f} '
          f'{ddjf:+7.2f} {gl:+8.2f} | {v}')

import csv as _csv
with open(f'{os.path.dirname(os.path.abspath(__file__))}/../../data/so_decision_table.csv',
          'w', newline='') as _f:
    _w = _csv.writer(_f)
    _w.writerow(['run', 'SO_sfc_bias', 'SO_sfc_rmse', 'SO_TOA_CRE_bias', 'dTropics',
                 'dSiberiaJJA', 'dSiberiaDJF', 'global_netTOA'])
    for r in rows:
        _w.writerow([r[0]] + [round(x, 3) for x in r[1:]])

print('\n  Ranked by |SO surface bias|.  "CLEAN" = tropics within +-0.50, and neither')
print('  Siberian season worse than its own threshold (JJA 0.246, DJF 0.604).')
print('  SO sfc / SO CRE are ABSOLUTE biases vs CERES; dTrop/dJJA/dDJF are vs control.')
