"""Why the Southern Ocean water column destabilised in PICAL_crunveg_kpp_gm1000.

The arm (KPP, uniform K_GM_max 1000) and the TKE-only line (cvmix_TKE, K_GM_max 2500 south of the
equator) share the ocean state of 2169-12-31. Within a few years the arm convects over about twice
the area south of 55S and loses 40 % of its winter sea ice. This follows the water column in the
region that newly convects, month by month from the branch point, in both runs.

Region NEW: open-ocean nodes (no cavity, bottom deeper than 2000 m) south of 55S whose September
mixed layer in 2177 is deeper than 1000 m in the arm and shallower than 500 m in the TKE-only line.
Region OLD (for reference): the same, but deeper than 1000 m in both (the Weddell convection the
line already has).

Per month, area-weighted over the region:
  S and T of the surface layer (0-50 m), of the winter-water layer (50-150 m) and of the warm deep
  water (200-500 m); the density step sigma0(200-300 m) - sigma0(0-50 m) and its split into a
  salinity part and a temperature part (linearised with alpha and beta of the mean state);
  Kv at 30 and 75 m (log mean); N2 max in the top 300 m; the GM bolus w at 200 m; surface
  freshwater flux fw and the part from sea ice fw_ice; net heat flux fh; ice concentration; MLD2.
The table gives the arm, the TKE-only line and arm minus TKE-only for selected months.

Usage: so_destabilisation.py     (reval environment, needs ~20 GB; writes data/clim/so_destabilisation*.{txt,csv})
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import pandas as pd
import xarray as xr
import gsw

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
RUNS = {'arm': f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_kpp_gm1000/outdata/fesom',
        'tke': f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082/outdata/fesom'}
Y0, Y1, YREF = 2170, 2179, 2177

md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lon, lat = md['lon'].values, md['lat'].values
zi = np.abs(md['nz'].values)
zm, dz = 0.5 * (zi[1:] + zi[:-1]), np.diff(zi)
area = md['nod_area'].values[0]
open_deep = (md['ulevels_nod2D'].values == 1) & (np.abs(md['zbar_n_bottom'].values) > 2000) & (lat < -55)


def sep_mld(run, y):
    with xr.open_dataset(f'{RUNS[run]}/MLD2.fesom.{y}.nc') as ds:
        return np.abs(ds['MLD2'].isel(time=8).values)


ma, mt = sep_mld('arm', YREF), sep_mld('tke', YREF)
REG = {'NEW': open_deep & (ma > 1000) & (mt < 500), 'OLD': open_deep & (ma > 1000) & (mt > 1000)}
LAY = {'0-50': (0, 50), '50-150': (50, 150), '200-300': (200, 300), '200-500': (200, 500)}


def lay_mean(x, a, b):
    """x: (nodes, nz) -> thickness-weighted layer mean per node."""
    m = (zm >= a) & (zm < b)
    return np.nansum(x[:, m] * dz[m], 1) / np.nansum(np.where(np.isfinite(x[:, m]), dz[m], 0.0), 1)


def at_depth(x, z, zax):
    return x[:, int(np.argmin(np.abs(zax - z)))]


rows = []
for run, d in RUNS.items():
    for y in range(Y0, Y1 + 1):
        ds = {v: xr.open_dataset(f'{d}/{v}.fesom.{y}.nc')[v] for v in ('temp', 'salt', 'Kv', 'N2', 'bolus_w', 'fw', 'fw_ice', 'fh', 'MLD2')}
        ai = xr.open_dataset(f'{d}/a_ice.fesom.{y}.nc')['a_ice']
        ai = ai.groupby('time.month').mean('time') if ai.sizes['time'] > 12 else ai
        for name, reg in REG.items():
            idx = np.where(reg)[0]
            w = area[idx] / area[idx].sum()
            for mo in range(12):
                T = ds['temp'].isel(time=mo).values[idx]
                S = ds['salt'].isel(time=mo).values[idx]
                T = np.where(S > 1, T, np.nan); S = np.where(S > 1, S, np.nan)
                r = dict(run=run, region=name, year=y, month=mo + 1)
                for k, (a, b) in LAY.items():
                    r[f'S{k}'] = float(np.nansum(lay_mean(S, a, b) * w))
                    r[f'T{k}'] = float(np.nansum(lay_mean(T, a, b) * w))
                s_top, s_bot = r['S0-50'], r['S200-300']
                t_top, t_bot = r['T0-50'], r['T200-300']
                r['dsig'] = float(gsw.sigma0(s_bot, t_bot) - gsw.sigma0(s_top, t_top))
                sm, tm = 0.5 * (s_top + s_bot), 0.5 * (t_top + t_bot)
                rho = 1000 + gsw.sigma0(sm, tm)
                r['dsig_S'] = float(rho * gsw.beta(sm, tm, 150) * (s_bot - s_top))
                r['dsig_T'] = float(-rho * gsw.alpha(sm, tm, 150) * (t_bot - t_top))
                kv = ds['Kv'].isel(time=mo).values[idx]
                for z in (30, 75):
                    r[f'logKv{z}'] = float(np.nansum(np.log10(np.maximum(at_depth(kv, z, zi), 1e-8)) * w))
                n2 = ds['N2'].isel(time=mo).values[idx][:, zi <= 300]
                r['N2max'] = float(np.nansum(np.nanmax(n2, 1) * w))
                r['bolus_w200'] = float(np.nansum(at_depth(ds['bolus_w'].isel(time=mo).values[idx], 200, zi) * w))
                for v in ('fw', 'fw_ice', 'fh'):
                    r[v] = float(np.nansum(ds[v].isel(time=mo).values[idx] * w))
                r['a_ice'] = float(np.nansum(ai.isel(month=mo).values[idx] * w)) if 'month' in ai.dims else float(np.nansum(ai.isel(time=mo).values[idx] * w))
                r['mld'] = float(np.nansum(np.abs(ds['MLD2'].isel(time=mo).values[idx]) * w))
                rows.append(r)
        print(run, y, flush=True)
D = pd.DataFrame(rows)
D.to_csv(f'{REPO}/data/clim/so_destabilisation_monthly.csv', index=False, float_format='%.6g')

SEC = 86400 * 30.0
COLS = [('S0-50', 'S 0-50 m', 1, '%8.3f'), ('S50-150', 'S 50-150 m', 1, '%8.3f'), ('S200-500', 'S 200-500 m', 1, '%8.3f'),
        ('T0-50', 'T 0-50 m', 1, '%8.2f'), ('T50-150', 'T 50-150 m', 1, '%8.2f'), ('T200-500', 'T 200-500 m', 1, '%8.2f'),
        ('dsig', 'sigma0 step', 1, '%8.3f'), ('dsig_S', '  salt part', 1, '%8.3f'), ('dsig_T', '  temp part', 1, '%8.3f'),
        ('logKv30', 'log10 Kv 30 m', 1, '%8.2f'), ('logKv75', 'log10 Kv 75 m', 1, '%8.2f'), ('N2max', 'N2 max 1e-5', 1e5, '%8.2f'),
        ('bolus_w200', 'bolus w 200m m/yr', 86400 * 365.0, '%8.1f'), ('fw', 'fw m/month', SEC, '%8.3f'), ('fw_ice', 'fw_ice m/month', SEC, '%8.3f'),
        ('fh', 'fh W/m2', 1, '%8.1f'), ('a_ice', 'ice conc', 1, '%8.2f'), ('mld', 'MLD2 m', 1, '%8.0f')]
WHEN = [(2170, 1), (2170, 3), (2170, 6), (2170, 9), (2171, 3), (2171, 9), (2172, 3), (2172, 9), (2173, 9), (2175, 9), (2177, 9), (2179, 9)]
out = [f"Regions from September {YREF}: NEW = MLD2 > 1000 m in the arm and < 500 m in the TKE-only line "
       f"({area[REG['NEW']].sum() / 1e12:.2f}e6 km2); OLD = > 1000 m in both ({area[REG['OLD']].sum() / 1e12:.2f}e6 km2). "
       "Open ocean south of 55S, bottom deeper than 2000 m.",
       'sigma0 step = sigma0(200-300 m) - sigma0(0-50 m) [kg/m3]; fw positive = freshwater into the ocean; fh positive = heat out of the ocean (FESOM sign as written).']
for name in REG:
    for lab, sel in (('arm (KPP, GM 1000)', 'arm'), ('TKE-only line', 'tke'), ('arm minus TKE-only', None)):
        out += ['', f'Region {name}: {lab}', f"{'quantity':18s}" + ''.join(f'{y}-{m:02d}'.rjust(9) for y, m in WHEN)]
        for key, nm, sc, fmt in COLS:
            line = f'{nm:18s}'
            for y, m in WHEN:
                q = D[(D.region == name) & (D.year == y) & (D.month == m)].set_index('run')[key] * sc
                val = q['arm'] - q['tke'] if sel is None else q[sel]
                line += ' ' + (fmt % val)
            out.append(line)
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/so_destabilisation.txt', 'w').write(txt + '\n')
