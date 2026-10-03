"""Why does the 1800 s FESOM step weaken Labrador convection?  Pathway comparison.

The first winter convects equally deep at 1200 and 1800 s (2-month tests perfE/G vs
perfF/H), so the step acts through preconditioning that builds over years.  Per run,
2121-2129 means of every pathway into the Labrador interior:
  ice drift      mean daily |u_ice| over a_ice >= 0.15 (central Arctic, Baffin, Labrador)
  ice budget     Nov-Apr dyngrice (import) and thdgrice (melt) in the interior
  surface forcing Nov-Mar fh (heat flux, W/m2, + into ocean) and fw (freshwater flux) in
                 the interior
  gyre doming    Nov-Mar curl of the total ocean surface stress (tx_sur, ty_sur) over the
                 Labrador basin, 1e-7 N/m3 (positive = cyclonic = upwelling, preconditions)
  preconditioning October temperature and salinity in the interior at 0-100, 100-500,
                 500-1500 m (native mesh)
Usage: python labrador_pathways.py
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4'
RUNS = [('PICAL_crunveg_ob1200', 'old 1200'), ('PICAL_crunveg_ob1800', 'old 1800'),
        ('PICAL_crunveg_ts1200', 'new 1200'), ('PICAL_crunveg_ihf0', 'new 1800')]
Y = range(2121, 2130)
INT = (56, 62, 300, 310); BASIN = (52, 64, 296, 315)
ICEB = {'>80N': (80, 90, 0, 360), 'Baffin': (66, 78, 280, 300), 'Labrador': (52, 66, 290, 316)}
NOVMAR = [10, 11, 0, 1, 2]; NOVAPR = [10, 11, 0, 1, 2, 3]
SEC = 86400 * 30.44


def gr(exp, v, y):
    with xr.open_dataset(f'{R}/{exp}/outdata/fesom/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
        a = np.squeeze(d[v].values).astype('f8'); lat = d['lat'].values; lon = d['lon'].values % 360
    return np.where(np.abs(a) > 1e30, np.nan, a), lat, lon


def monthly(a):
    if a.shape[0] == 12:
        return a
    ml = [31, 29 if a.shape[0] == 366 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    e = np.cumsum([0] + ml); return np.stack([np.nanmean(a[e[i]:e[i + 1]], 0) for i in range(12)])


def box(LAT, LON, b):
    return (LAT >= b[0]) & (LAT <= b[1]) & (LON >= b[2]) & (LON <= b[3])


res = {}
for exp, lab in RUNS:
    r = {}
    for y in Y:
        a, lat, lon = gr(exp, 'a_ice', y)
        LAT, LON = np.meshgrid(lat, lon, indexing='ij'); W = np.cos(np.deg2rad(LAT))
        u, _, _ = gr(exp, 'uice', y); v, _, _ = gr(exp, 'vice', y)
        if u.shape[0] == a.shape[0]:
            s = np.hypot(u, v)
            for n, b in ICEB.items():
                m = (a >= 0.15) & box(LAT, LON, b)[None]
                r.setdefault(f'drift {n} cm/s', []).append(100 * np.nansum(np.where(m, s, 0) * W) / np.sum(m * W))
        am = monthly(a); kI = box(LAT, LON, INT); kB = box(LAT, LON, BASIN)
        f = lambda x, k, ms: np.mean([np.nansum((np.nan_to_num(x[m]) * W)[k]) / W[k].sum() for m in ms])
        for v_, nm, ms, sc in (('dyngrice', 'ice import cm', NOVAPR, SEC * 100 * 6), ('thdgrice', 'ice thermo cm', NOVAPR, SEC * 100 * 6),
                               ('fh', 'fh W/m2', NOVMAR, 1), ('fw', 'fw mm/day', NOVMAR, 86400e3)):
            try:
                x, _, _ = gr(exp, v_, y); r.setdefault(nm, []).append(f(monthly(x), kI, ms) * sc)
            except Exception:
                pass
        try:
            tx, _, _ = gr(exp, 'tx_sur', y); ty, _, _ = gr(exp, 'ty_sur', y)
            tx, ty = monthly(tx), monthly(ty)
            dy = np.deg2rad(lat[1] - lat[0]) * 6.371e6; dx = np.deg2rad(lon[1] - lon[0]) * 6.371e6 * np.cos(np.deg2rad(LAT))
            curl = np.gradient(ty, axis=2) / dx - np.gradient(tx, axis=1) / dy
            r.setdefault('stress curl 1e-7 N/m3', []).append(f(curl, kB, NOVMAR) * 1e7)
        except Exception:
            pass
    # October profiles in the interior, native mesh
    for v_ in ('salt', 'temp'):
        prof = []
        for y in Y:
            with xr.open_dataset(f'{R}/{exp}/outdata/fesom/{v_}.fesom.{y}.nc', decode_times=False) as d:
                x = d[v_].isel(time=9).values; dz = [c for c in d[v_].dims if c.startswith('nz')][0]
                dep = np.abs(d[dz].values); lo = ((d['lon'].values + 180) % 360) - 180; la = d['lat'].values
            if x.shape[0] == len(dep): x = x.T
            k = (la >= 56) & (la < 62) & (lo >= -60) & (lo < -50)
            x = x[k]; x = np.where(np.abs(x) > 1e30, np.nan, x)
            prof.append([np.nanmean(x[:, (dep >= z0) & (dep < z1)]) for z0, z1 in ((0, 100), (100, 500), (500, 1500))])
        for i, z in enumerate(('0-100', '100-500', '500-1500')):
            r[f'Oct {v_} {z} m'] = [p[i] for p in prof]
    res[lab] = r
keys = list(res['old 1200'].keys())
print(f'{"2121-2129 mean (sd)":28s}' + ''.join(f'{l:>18s}' for _, l in RUNS))
for k in keys:
    print(f'{k:28s}' + ''.join(f'{np.mean(res[l][k]):10.3f} ({np.std(res[l][k], ddof=1):5.3f})' if k in res[l] else f'{"--":>18s}' for _, l in RUNS))
