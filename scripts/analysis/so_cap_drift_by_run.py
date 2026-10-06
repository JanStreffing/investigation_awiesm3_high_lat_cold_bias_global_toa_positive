"""Southern Ocean cap and warm deep water in every run with 3D output on albedo: state and drift.

Fixed open-ocean region: south of 55S, bottom deeper than 2000 m, no cavity (15e6 km2 scale), and the
two sub-regions NEW / OLD of so_destabilisation.py. For every second year of each run, September:
temperature and salinity of 0-50 m, 50-150 m (winter water) and 200-500 m (warm deep water), the
density step sigma0(200-300 m) - sigma0(0-50 m), and the annual-mean temperature of 200-500 m and
500-1500 m. Per run: mean of the first and of the last three sampled years and the linear trend per
decade. PHC3 September values of the same regions are in data/clim/so_cap_vs_phc.txt.

Usage: so_cap_drift_by_run.py    (reval environment; writes data/clim/so_cap_drift_by_run.{txt,csv})
"""
import os
import numpy as np
import pandas as pd
import xarray as xr
import gsw

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
B = f'{P}/runtime/awiesm3-v3.4'
RUNS = [('ccnice 1200s GM2500 TKE+IDEMIX', 'PICAL_ccnice', 2101, 2129),
        ('crunveg 1200s GM2500 TKE+IDEMIX', 'PICAL_crunveg', 2101, 2119),
        ('gm2500 1200s (ob1200)', 'PICAL_crunveg_ob1200', 2121, 2129),
        ('gm1500 1200s', 'PICAL_crunveg_gm1500', 2121, 2129),
        ('gmramp 1200s', 'PICAL_crunveg_gmramp', 2121, 2129),
        ('ob1800 GM2500', 'PICAL_crunveg_ob1800', 2121, 2139),
        ('gmferr1800 Ferreira', 'PICAL_crunveg_gmferr1800', 2121, 2139),
        ('gmrd1800 Rossby radius', 'PICAL_crunveg_gmrd1800', 2121, 2139),
        ('gmhemi1800 TKE+IDEMIX', 'PICAL_crunveg_gmhemi1800', 2121, 2139),
        ('gmhemi_tke TKE', 'PICAL_crunveg_gmhemi_tke', 2141, 2159),
        ('tke_albsn082 TKE', 'PICAL_crunveg_tke_albsn082', 2161, 2209),
        ('idemix_albsn082 TKE+IDEMIX', 'PICAL_crunveg_idemix_albsn082', 2171, 2209),
        ('kpp_gm1000 KPP GM1000', 'PICAL_crunveg_kpp_gm1000', 2171, 2189)]
md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lat = md['lat'].values
zi = np.abs(md['nz'].values); zm, dz = 0.5 * (zi[1:] + zi[:-1]), np.diff(zi)
area = md['nod_area'].values[0]
od = (md['ulevels_nod2D'].values == 1) & (np.abs(md['zbar_n_bottom'].values) > 2000) & (lat < -55)
mld = lambda e: np.abs(xr.open_dataset(f'{B}/{e}/outdata/fesom/MLD2.fesom.2177.nc')['MLD2'].isel(time=8).values)
ma, mt = mld('PICAL_crunveg_kpp_gm1000'), mld('PICAL_crunveg_tke_albsn082')
REG = {'SO55': od, 'NEW': od & (ma > 1000) & (mt < 500), 'OLD': od & (ma > 1000) & (mt > 1000)}
IDX = {k: np.where(v)[0] for k, v in REG.items()}
W = {k: area[i] / area[i].sum() for k, i in IDX.items()}
LAY = [('0-50', 0, 50), ('50-150', 50, 150), ('200-300', 200, 300), ('200-500', 200, 500), ('500-1500', 500, 1500)]


def lay(x, a, b):
    m = (zm >= a) & (zm < b)
    return np.nansum(x[:, m] * dz[m], 1) / np.nansum(np.where(np.isfinite(x[:, m]), dz[m], 0), 1)


rows = []
for label, exp, y0, y1 in RUNS:
    d = f'{B}/{exp}/outdata/fesom'
    for y in range(y0, y1 + 1, 2):
        try:
            Tt = xr.open_dataset(f'{d}/temp.fesom.{y}.nc')['temp']; St = xr.open_dataset(f'{d}/salt.fesom.{y}.nc')['salt']
            T9, S9, Ta = Tt.isel(time=8).values, St.isel(time=8).values, Tt.mean('time').values
        except Exception as e:
            print(label, y, 'skipped', str(e)[:60], flush=True); continue
        ok = S9 > 1
        T9 = np.where(ok, T9, np.nan); S9 = np.where(ok, S9, np.nan); Ta = np.where(ok, Ta, np.nan)
        for k, idx in IDX.items():
            r = dict(run=label, region=k, year=y)
            for nm, a, b in LAY:
                r[f'T{nm}'] = float(np.nansum(lay(T9[idx], a, b) * W[k])); r[f'S{nm}'] = float(np.nansum(lay(S9[idx], a, b) * W[k]))
                r[f'Tann{nm}'] = float(np.nansum(lay(Ta[idx], a, b) * W[k]))
            r['step'] = float(gsw.sigma0(r['S200-300'], r['T200-300']) - gsw.sigma0(r['S0-50'], r['T0-50']))
            rows.append(r)
        print(label, y, flush=True)
D = pd.DataFrame(rows)
D.to_csv(f'{REPO}/data/clim/so_cap_drift_by_run.csv', index=False, float_format='%.6g')


def tr(y, x):
    return 10 * np.polyfit(y, x, 1)[0] if len(y) > 2 else np.nan


out = ['September state and drift of the Southern Ocean cap by run (every second year). first/last = mean of the first/last three',
       'sampled years; /dec = linear trend per decade. step = sigma0(200-300 m) - sigma0(0-50 m). Tann = annual mean.',
       'PHC3 September: region NEW step 0.320, S 0-50 34.04, T 50-150 -1.12, T 200-500 0.97; region OLD step 0.238, S 0-50 34.22, T 50-150 -1.10, T 200-500 0.58.']
for k in REG:
    out += ['', f'Region {k} ({area[IDX[k]].sum() / 1e12:.2f}e6 km2)',
            f"{'run':32s}{'years':>11s}{'step first':>11s}{'last':>7s}{'/dec':>8s}{'S0-50 last':>11s}{'T50-150 last':>13s}{'/dec':>7s}"
            f"{'Tann200-500 first':>18s}{'last':>7s}{'/dec':>8s}{'Tann500-1500 /dec':>18s}"]
    for label, exp, y0, y1 in RUNS:
        q = D[(D.run == label) & (D.region == k)].sort_values('year')
        if len(q) < 3:
            continue
        y = q.year.values.astype(float)
        f3, l3 = q.iloc[:3], q.iloc[-3:]
        out.append(f"{label:32s}{f'{int(y[0])}-{int(y[-1])}':>11s}{f3.step.mean():11.3f}{l3.step.mean():7.3f}{tr(y, q.step.values):+8.3f}"
                   f"{l3['S0-50'].mean():11.3f}{l3['T50-150'].mean():13.2f}{tr(y, q['T50-150'].values):+7.2f}"
                   f"{f3['Tann200-500'].mean():18.3f}{l3['Tann200-500'].mean():7.3f}{tr(y, q['Tann200-500'].values):+8.3f}{tr(y, q['Tann500-1500'].values):+18.3f}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/so_cap_drift_by_run.txt', 'w').write(txt + '\n')
