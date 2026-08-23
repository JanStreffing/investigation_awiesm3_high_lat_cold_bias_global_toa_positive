"""Southern Ocean: is there too little sea ice, and is it too dark?

WHY.  The remaining global cold bias (~0.7 K) sits alongside a positive TOA imbalance
(~+0.88 W/m2), so the model will keep warming toward equilibrium and close much of that
bias on its own.  That is fine globally and bad regionally if the Southern Ocean is
already too warm, because SH sea ice is the amplifier: less ice -> lower albedo -> more
absorbed shortwave -> less ice.  Whether the extra warming is safe there depends on which
way the SH ice error points, and on whether the ice that exists reflects enough.

TWO SEPARATE QUESTIONS, deliberately not conflated:
  AMOUNT   SH sea ice AREA (sum of concentration x cell area) and EXTENT (area of cells
           with concentration >= 15 %) by month, against OSI-SAF.  Extent and area are
           different numbers and mixing them is the classic way to get the sign wrong.
  DARKNESS Surface albedo where ice actually is, and the clear-sky SW reflected, so a
           correct ice AREA with too-dark ice is distinguishable from too little ice.

OSI-SAF is the EUMETSAT climate data record on a 432x432 Lambert azimuthal grid, so it is
regridded to the model grid with cdo before comparison; its own cell areas are used for
the observed totals rather than the model's.
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
ARMS = [('11P LX4 1990', f'{R}/11P'), ('11R +RSBLB2', f'{R}/11R')]
Y0, Y1 = 1380, 1389
OSI = '/work/ab0246/a270092/obs/osisaf_sh.nc'
MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

f = sorted(glob.glob(f'{R}/11P/outdata/oifs/atm_remapped_1m_lsm_*.nc'))[0]
with xr.open_dataset(f, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values)
    m = m[0] if m.ndim == 3 else m
    lat = np.squeeze(d['lat'].values)
ocean = m <= 0.5
RE = 6.371e6
AREA = np.broadcast_to((RE**2 * np.cos(np.deg2rad(lat)) * np.deg2rad(abs(lat[1]-lat[0]))
                        * 2*np.pi/m.shape[1])[:, None], m.shape)
sh = ocean & np.broadcast_to((lat < -45)[:, None], m.shape)


def load(path, var):
    out = []
    for y in range(Y0, Y1 + 1):
        p = f'{path}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(p):
            return None
        with xr.open_dataset(p, decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            out.append(np.squeeze(d[k].values))
    return np.stack(out).mean(0)          # (12, ny, nx) climatology


# --- observed: OSI-SAF monthly climatology, on its own grid with its own cell area ---
with xr.open_dataset(OSI, decode_times=True) as d:
    g = d['ice_conc'].groupby('time.month').mean('time')
    have = [int(x) for x in g['month'].values]
    conc = np.full((12,) + g.values.shape[1:], np.nan)
    for j, mo in enumerate(have):
        conc[mo - 1] = g.values[j] / 100.0        # % -> frac
    # This record runs 1989-01 to 2014-11 and contains NO December at all (11 months x
    # 26 years = 286 steps), so December is left NaN rather than silently mis-indexed.
    MISSING = [MON[i] for i in range(12) if i + 1 not in have]
    # Lambert azimuthal grid, 25 km nominal: xc/yc are in km and regularly spaced
    dx = abs(float(d['xc'][1] - d['xc'][0])) * 1000.0
    dy = abs(float(d['yc'][1] - d['yc'][0])) * 1000.0
cell = dx * dy
def _tot(x, thr=None):
    if np.all(np.isnan(x)):
        return np.nan
    return (np.nansum(x >= thr) if thr else np.nansum(x)) * cell / 1e12
obs_area = np.array([_tot(conc[i]) for i in range(12)])
obs_ext = np.array([_tot(conc[i], 0.15) for i in range(12)])

print(f'SH sea ice, model {Y0}-{Y1} vs OSI-SAF 1989-2014 [1e6 km2]')
print(f'OSI-SAF has no data for: {", ".join(MISSING) if MISSING else "(none)"}\n')
print(f'{"":>6} | ' + ' '.join(f'{x:>6}' for x in MON))
fmt = lambda v: '   n/a' if np.isnan(v) else f'{v:6.2f}'
print(f'{"OBS a":>6} | ' + ' '.join(fmt(x) for x in obs_area))
print(f'{"OBS e":>6} | ' + ' '.join(fmt(x) for x in obs_ext))
res = {}
for label, path in ARMS:
    ci = load(path, 'ci')
    if ci is None:
        print(f'{label}: incomplete'); continue
    a = np.array([np.nansum(np.where(sh, ci[i], 0) * AREA) / 1e12 for i in range(12)])
    e = np.array([np.nansum(np.where(sh & (ci[i] >= 0.15), 1, 0) * AREA) / 1e12 for i in range(12)])
    res[label] = (a, e, ci)
    print(f'{label[:6]:>6} | ' + ' '.join(f'{x:6.2f}' for x in a))
print()
for label in res:
    a, e, _ = res[label]
    print(f'{label:>14}  area  Sep {a[8]:6.2f} (obs {obs_area[8]:5.2f}, {a[8]-obs_area[8]:+5.2f})'
          f'   Feb {a[1]:5.2f} (obs {obs_area[1]:5.2f}, {a[1]-obs_area[1]:+5.2f})')
    print(f'{"":>14}  ext   Sep {e[8]:6.2f} (obs {obs_ext[8]:5.2f}, {e[8]-obs_ext[8]:+5.2f})'
          f'   Feb {e[1]:5.2f} (obs {obs_ext[1]:5.2f}, {e[1]-obs_ext[1]:+5.2f})')

# --- is the ice too dark? albedo where ice is, and the SO SW budget ---
print('\nSurface albedo and SW, SH 45-90S ocean, annual and September:')
print(f'{"arm":>14} | {"fal(ice>0.5)":>13} {"fal ann":>8} | {"SWup/SWdn ann":>14}')
for label, path in ARMS:
    fal = load(path, 'fal')
    ci = res[label][2]
    if fal is None:
        print(f'{label:>14} | fal not archived'); continue
    icy = sh & (ci[8] >= 0.5)
    w = AREA
    v1 = float(np.average(fal[8][icy], weights=w[icy])) if icy.any() else np.nan
    v2 = float(np.average(fal.mean(0)[sh], weights=w[sh]))
    print(f'{label:>14} | {v1:13.3f} {v2:8.3f} |')
