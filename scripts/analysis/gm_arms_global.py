"""What do the GM tests do outside the Labrador Sea? Annual global diagnostics for
PICAL_crunveg_gm1500 / gmramp against their reference ob1200 (same build, step, start;
GM 2500 without ramp), plus the ccnice continuation for context.

Per run and year (2120 is the OpenIFS cold start and is skipped):
  T2m global / NH / SH / 60-90N / 60-90S / land / Siberia (55-75N 60-180E land; annual, DJF, JJA)
  TOA net (tsr+ttr), net surface flux into the surface (ssr+str+slhf+sshf - Lf*snowfall), from the remapped
    monthly OIFS fields (accumulated over the 1 h output step, so /3600 -> W m-2)
  FESOM net heat flux into the ocean (-fh, area-weighted)
  ocean heat content by depth band (0-700, 700-2000, 2000+ m) and by region (SO < 40S,
    SPNA 45-66N 80W-10E, NORDARC > 66N, REST), as the rate of change over the period
  sea-ice volume NH March / September, SH September / February; extent NH March / Sep, SH Sep / Feb
Usage: python gm_arms_global.py <y0> <y1> tag=fesomdir:oifsdir [...]   (first tag is the reference;
       call PICAL_crunveg_ob1200 "gm2500")
Output: data/clim/gm_arms_global_<y0>-<y1>.csv (annual) and a table vs the first run.
"""
import os, sys, numpy as np, pandas as pd, xarray as xr
from netCDF4 import Dataset
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
y0, y1 = int(sys.argv[1]), int(sys.argv[2])
RUNS = [(a.split('=')[0], *a.split('=')[1].split(':')) for a in sys.argv[3:]]
LSM = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4/PICAL_ccnice/outdata/oifs/atm_remapped_1m_lsm_2110-2110.nc'
DIAG = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/input/fesom2/core3/fesom.mesh.diag.nc'
with xr.open_dataset(LSM) as d:
    lsm = d.lsm.isel(time_counter=0).load() if 'time_counter' in d.lsm.dims else d.lsm.load()
with Dataset(DIAG) as nc:
    area = np.asarray(nc.variables['nod_area'][:], 'f8'); zb = np.abs(np.asarray(nc.variables['nz'][:], 'f8'))
    mlon = np.asarray(nc.variables['lon'][:], 'f8'); mlat = np.asarray(nc.variables['lat'][:], 'f8')
if np.abs(mlon).max() < 7: mlon, mlat = np.rad2deg(mlon), np.rad2deg(mlat)
mlon = ((mlon + 180) % 360) - 180
zmid = 0.5 * (zb[:-1] + zb[1:])
REG = {'SO': mlat < -40, 'SPNA': (mlat >= 45) & (mlat < 66) & (mlon >= -80) & (mlon <= 10), 'NORDARC': mlat >= 66}
REG['REST'] = ~(REG['SO'] | REG['SPNA'] | REG['NORDARC'])
BAND = {'0-700': zmid < 700, '700-2000': (zmid >= 700) & (zmid < 2000), '2000+': zmid >= 2000}
AE, YR = 5.1e14, 365.25 * 86400

def oifs(d, v, y):
    f = f'{d}/atm_remapped_1m_{v}_{y}-{y}.nc'
    if not os.path.exists(f): return None
    with xr.open_dataset(f) as ds: return ds[v].load()

def wm(x, w, m=None):
    w = w if m is None else w.where(m, 0.0)
    return float((x * w).sum(('lat', 'lon')) / w.sum(('lat', 'lon')))

def fes(d, v, y):
    f = f'{d}/{v}.fesom.{y}.nc'
    if not os.path.exists(f): return None
    with Dataset(f) as nc:
        x = nc.variables[v][:]
        return np.asarray(x.filled(np.nan) if hasattr(x, 'filled') else x, 'f8')

rows = []
for tag, fd, od in RUNS:
    h = fes(fd, 'hnode', y0 - 1 if os.path.exists(f'{fd}/hnode.fesom.{y0-1}.nc') else y0)
    h = np.nanmean(h, 0); h = h if h.shape[0] == 47 else h.T
    vol = np.where(np.isfinite(h), h, 0) * area[:47]
    for y in range(y0 - 1, y1 + 1):          # y0-1 only to anchor the OHC change
        r = {'run': tag, 'year': y}
        t = fes(fd, 'temp', y)
        if t is not None:
            t = np.nanmean(np.where(np.abs(t) > 1e5, np.nan, t), 0); t = t if t.shape[0] == 47 else t.T
            hc = 4.1e6 * np.nan_to_num(t) * vol
            for b, k in BAND.items(): r[f'ohc_{b}'] = hc[k].sum()
            for g, m in REG.items(): r[f'ohc_{g}'] = hc[:, m].sum()
        if y < y0: rows.append(r); continue
        t2 = oifs(od, '2t', y)
        if t2 is None: print(tag, y, 'no 2t', flush=True); rows.append(r); continue
        w = np.cos(np.deg2rad(t2.lat)) * xr.ones_like(t2.lon)
        lon = t2.lon % 360
        ann = t2.mean('time_counter') - 273.15
        land = xr.DataArray(np.squeeze(lsm.values) > 0.5, coords=ann.coords, dims=ann.dims)
        sib = land & (t2.lat >= 55) & (t2.lat <= 75) & (lon >= 60) & (lon <= 180)
        djf = t2.isel(time_counter=[0, 1, 11]).mean('time_counter') - 273.15
        jja = t2.isel(time_counter=[5, 6, 7]).mean('time_counter') - 273.15
        r.update(t2m_glob=wm(ann, w), t2m_nh=wm(ann, w, t2.lat >= 0), t2m_sh=wm(ann, w, t2.lat < 0),
                 t2m_6090n=wm(ann, w, t2.lat >= 60), t2m_6090s=wm(ann, w, t2.lat <= -60),
                 t2m_land=wm(ann, w, land), t2m_sib=wm(ann, w, sib), t2m_sib_djf=wm(djf, w, sib),
                 t2m_sib_jja=wm(jja, w, sib))
        F = {v: oifs(od, v, y) for v in ('tsr', 'ttr', 'ssr', 'str', 'slhf', 'sshf', 'sf')}
        if all(F[v] is not None for v in F):
            g = {v: wm(F[v].mean('time_counter'), w) / 3600.0 for v in F}
            r['toa_net'] = g['tsr'] + g['ttr']
            # the four radiative/turbulent terms omit the latent heat of fusion that falling
            # snow takes from the surface (sf in m water per step): about 0.83 W m-2 globally
            r['sfc_net'] = g['ssr'] + g['str'] + g['slhf'] + g['sshf'] - 3.337e5 * 1000.0 * g['sf']
        fh = fes(fd, 'fh', y)
        if fh is not None:
            fh = np.nanmean(np.where(np.abs(fh) > 1e5, np.nan, fh), 0); ok = np.isfinite(fh)
            r['ocean_in_Wm2_earth'] = -np.sum(np.where(ok, fh, 0) * area[0]) / AE
        for v, m in (('sivoln', 2), ('sivoln', 8), ('sivols', 1), ('sivols', 8),
                     ('siextentn', 2), ('siextentn', 8), ('siextents', 1), ('siextents', 8)):
            x = fes(fd, v, y)
            if x is not None and x.size == 12: r[f'{v}_m{m+1:02d}'] = float(np.ravel(x)[m])
        rows.append(r); print(tag, y, 'ok', flush=True)
df = pd.DataFrame(rows).sort_values(['run', 'year'])
for c in [c for c in df.columns if c.startswith('ohc_')]:      # J -> W m-2 of Earth surface
    df[c.replace('ohc_', 'heat_')] = df.groupby('run')[c].diff() / YR / AE
out = f'data/clim/gm_arms_global_{y0}-{y1}.csv'; df.to_csv(out, index=False)
d = df[df.year >= y0]
ref = RUNS[0][0]
cols = ['t2m_glob', 't2m_nh', 't2m_sh', 't2m_6090n', 't2m_6090s', 't2m_land', 't2m_sib', 't2m_sib_djf', 't2m_sib_jja',
        'toa_net', 'sfc_net', 'ocean_in_Wm2_earth', 'heat_0-700', 'heat_700-2000', 'heat_2000+',
        'heat_SO', 'heat_SPNA', 'heat_NORDARC', 'heat_REST',
        'sivoln_m03', 'sivoln_m09', 'sivols_m02', 'sivols_m09', 'siextentn_m03', 'siextentn_m09', 'siextents_m02', 'siextents_m09']
cols = [c for c in cols if c in d.columns]
M = d.groupby('run')[cols].mean().T; S = d.groupby('run')[cols].std().T; N = d.groupby('run')[cols].count().T
print(f'\n{y0}-{y1} means; differences vs {ref}; "*" = beyond the paired 2 sigma/sqrt(N) of the two runs')
hdr = f'{"":22s}' + ''.join(f'{t:>14s}' for t, *_ in RUNS) + ''.join(f'{t + "-" + ref:>22s}' for t, *_ in RUNS[1:])
print(hdr)
for c in cols:
    line = f'{c:22s}' + ''.join(f'{M.loc[c, t]:14.3f}' if t in M.columns else f'{"":14s}' for t, *_ in RUNS)
    for t, *_ in RUNS[1:]:
        if t not in M.columns: continue
        dd = M.loc[c, t] - M.loc[c, ref]
        thr = 2 * np.sqrt(S.loc[c, t]**2 / N.loc[c, t] + S.loc[c, ref]**2 / N.loc[c, ref])
        line += f'{dd:21.3f}{"*" if abs(dd) > thr else " "}'
    print(line)
