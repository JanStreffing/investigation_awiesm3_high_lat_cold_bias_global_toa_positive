"""Annual series of the v3.5.0-cc CMIP7 spin-up (fresh start from the PHC3 ocean-only restart, REcoM
cold, concentration-driven PI CO2, icebergs from the calving budget; monthly output only).

Per year, from the monthly files:
  atmosphere (remapped 400x192)  T2m global / 60-90N / 60-90S / land / ocean, net TOA and its two terms,
                                 global CO2 fluxes seen by OpenIFS (land co2nbf, ocean co2of), PgC/yr
  sea ice (native mesh)          area and volume, NH March / September, SH February / September
  ocean surface                  SST, SSS, air-sea CO2 flux (PgC/yr, positive into the ocean), net primary
                                 production (PgC/yr), iceberg melt (Gt/yr), September area of MLD2 > 1000 m
                                 south of 50S, March maximum MLD2 north of 45N in the Atlantic sector
  ocean interior                 volume means of temp, salt, DIC, Alk, O2, DIN in six depth bands and in total
  land (LPJ-GUESS, 1 degree)     total land carbon (PgC) and NEE (PgC/yr)

Usage: spinup_v350_annual_series.py <exp root> <y0> <y1> <out csv> [nproc]
(module analysis-toolbox/python-04.2026: its HDF5 reads the Blosc-compressed XIOS files)
"""
import os
import sys
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import glob
import numpy as np
import pandas as pd
import xarray as xr
from multiprocessing import Pool

ROOT, Y0, Y1, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
NP = int(sys.argv[5]) if len(sys.argv) > 5 else 8
MESH = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/reval_obs/mesh/core3/fesom.mesh.diag.nc'
BANDS = [(0, 100), (100, 300), (300, 700), (700, 1500), (1500, 2500), (2500, 6500)]
RE, YR = 6.371e6, 365.25 * 86400
TR3D = ('temp', 'salt', 'DIC', 'Alk', 'O2', 'DIN')

md = xr.open_dataset(MESH)
lat, lon = md['lat'].values, md['lon'].values
a3 = md['nod_area'].values                       # (nz, nod2)
area = a3[0]
zi = np.abs(md['nz'].values); zm, dz = 0.5 * (zi[1:] + zi[:-1]), np.diff(zi)
cav = md['ulevels_nod2D'].values > 1
vol = (a3[:len(dz)] * dz[:, None]).T             # (nod2, nz-1)
natl = (lat > 45) & (lon > -70) & (lon < 20) & ~cav


def f2(var, y):
    fn = f'{ROOT}/outdata/fesom/{var}.fesom.{y}.nc'
    if not os.path.isfile(fn):
        return None
    with xr.open_dataset(fn) as ds:
        return ds[var].values.astype(np.float64)


def oifs(var, y):
    fn = f'{ROOT}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
    if not os.path.isfile(fn):
        return None
    with xr.open_dataset(fn) as ds:
        return ds[var].load()


def gmean(x, w, m=None):
    w2 = w if m is None else w * m
    return float((x * w2).sum() / w2.sum())


def year(y):
    r = {'year': y}
    # atmosphere
    t = oifs('2t', y)
    if t is not None:
        w = np.cos(np.deg2rad(t['lat'])) + 0 * t.isel(time_counter=0)
        ta = t.mean('time_counter') - 273.15
        lsm = oifs('lsm', y).isel(time_counter=0)
        la = ta['lat']
        r.update(t2m_glob=gmean(ta, w), t2m_6090n=gmean(ta, w, la > 60), t2m_6090s=gmean(ta, w, la < -60),
                 t2m_land=gmean(ta, w, lsm > 0.5), t2m_ocean=gmean(ta, w, lsm <= 0.5),
                 t2m_6090n_djf=gmean(t.isel(time_counter=[0, 1, 11]).mean('time_counter') - 273.15, w, la > 60),
                 t2m_6090n_jja=gmean(t.isel(time_counter=[5, 6, 7]).mean('time_counter') - 273.15, w, la > 60))
        for v in ('tsr', 'ttr', 'ssr', 'str', 'sshf', 'slhf', 'sf'):
            x = oifs(v, y)
            if x is not None:
                r[v] = gmean(x.mean('time_counter'), w) / 3600.0
        if 'tsr' in r and 'ttr' in r:
            r['toa_net'] = r['tsr'] + r['ttr']
        if all(k in r for k in ('ssr', 'str', 'sshf', 'slhf', 'sf')):
            r['sfc_net'] = r['ssr'] + r['str'] + r['sshf'] + r['slhf'] - r['sf'] * 1000.0 * 3.3355e5   # sf: m w.e. per s after /3600
        aglob = 4 * np.pi * RE ** 2
        for v in ('co2nbf', 'co2of', 'co2apf'):
            x = oifs(v, y)
            if x is not None:
                r[v + '_PgC'] = gmean(x.mean('time_counter'), w) * aglob * YR * 12.0 / 44.0 / 1e12
    # sea ice
    a, m = f2('a_ice', y), f2('m_ice', y)
    if a is not None:
        for nm, sel, mo in (('nh_mar', lat > 40, 2), ('nh_sep', lat > 40, 8), ('sh_feb', lat < -50, 1), ('sh_sep', lat < -50, 8)):
            r['sia_' + nm] = np.nansum((a[mo] * area)[sel]) / 1e12
            r['siv_' + nm] = np.nansum((m[mo] * area)[sel]) / 1e12
    # surface ocean
    oc = ~cav
    for v in ('sst', 'sss'):
        x = f2(v, y)
        if x is not None:
            r[v] = float(np.nansum(x.mean(0)[oc] * area[oc]) / area[oc].sum())
    x = f2('CO2f', y)
    if x is not None:
        r['co2f_PgC'] = float(np.nansum(x.mean(0) * area)) * 365.25 * 12.011e-3 / 1e15
        r['co2f_PgC_s30'] = float(np.nansum((x.mean(0) * area)[lat < -30])) * 365.25 * 12.011e-3 / 1e15
    npp = [f2(v, y) for v in ('NPPn', 'NPPd')]
    if npp[0] is not None:
        r['npp_PgC'] = float(sum(np.nansum(p.mean(0) * area) for p in npp if p is not None)) * 365.25 * 12.011e-3 / 1e15
    ib = [f2(v, y) for v in ('ibfwb', 'ibfwl', 'ibfwe', 'ibfwbv')]
    if ib[0] is not None:
        r['icb_melt_Gt'] = float(sum(np.nansum(p.mean(0) * area) for p in ib)) * 1000.0 * YR / 1e12
    x = f2('MLD2', y)
    if x is not None:
        x = np.abs(x)
        r['so_mld1000_sep_Mkm2'] = float(area[(lat < -50) & oc & (x[8] > 1000)].sum()) / 1e12
        r['natl_mld_mar_max'] = float(np.nanmax(x[2][natl]))
        r['natl_mld1000_mar_Mkm2'] = float(area[natl & (x[2] > 1000)].sum()) / 1e12
    # interior
    S = f2('salt', y)
    if S is not None:
        wet = (S.mean(0) > 1.0)
        V = np.where(wet, vol, 0.0)
        for v in TR3D:
            x = S if v == 'salt' else f2(v, y)
            if x is None:
                continue
            xm = np.where(wet, x.mean(0), 0.0)
            r[f'{v}_tot'] = float((xm * V).sum() / V.sum())
            for (za, zb) in BANDS:
                k = (zm >= za) & (zm < zb)
                r[f'{v}_{za}-{zb}'] = float((xm[:, k] * V[:, k]).sum() / V[:, k].sum())
        r['ocean_vol'] = float(V.sum())
    return r


def land_fraction_1deg(lon_deg, lat_deg):
    """Fraction of each 1 degree cell (180 x 360, from -89.5 and -179.5) covered by the cells of the given
    points of the octahedral grid. The conservative remap of LPJ-GUESS gives a 1 degree cell the mean over
    its model-land part and does not keep how large that part is; without it totals come out 12 % high."""
    lat_deg, lon_deg = np.asarray(lat_deg, float), ((np.asarray(lon_deg, float) + 180.0) % 360.0) - 180.0
    rows_ = np.unique(np.round(lat_deg, 4))
    n = int(round(90.0 / np.min(np.diff(rows_))))
    xg, wg = np.polynomial.legendre.leggauss(2 * n)
    j = np.abs(np.rad2deg(np.arcsin(xg))[None, :] - lat_deg[:, None]).argmin(1)
    half = 180.0 / (20 + 4 * np.minimum(j, 2 * n - 1 - j))
    edge = np.concatenate([[-1.0], np.cumsum(wg) - 1.0]); edge[-1] = 1.0
    s0, s1, w, e = edge[j], edge[j + 1], lon_deg - half, lon_deg + half
    sin_e = np.sin(np.deg2rad(np.arange(-90.0, 90.5, 1.0)))
    acc = np.zeros((180, 360))
    for p in range(len(s0)):
        j0 = max(np.searchsorted(sin_e, s0[p], 'right') - 1, 0)
        j1 = min(np.searchsorted(sin_e, s1[p], 'left'), 180)
        dsin = np.clip(np.minimum(sin_e[j0 + 1:j1 + 1], s1[p]) - np.maximum(sin_e[j0:j1], s0[p]), 0, None)
        for i in range(int(np.floor(w[p] + 180.0)), int(np.ceil(e[p] + 180.0))):
            dlon = min(i + 1 - 180.0, e[p]) - max(i - 180.0, w[p])
            if dlon > 0:
                acc[j0:j1, i % 360] += dsin * dlon
    return np.clip(acc / np.diff(sin_e)[:, None], 0.0, 1.0)


def lpjg():
    rows = {}
    # the model's land points, from any leg that kept the native table
    frac = None
    for fn in sorted(glob.glob(f'{ROOT}/outdata/lpj_guess/*-*/cpool.out')):
        t = pd.read_csv(fn, sep=r'\s+')
        t = t[t['Year'] == t['Year'].iloc[0]]
        frac = land_fraction_1deg(t['Lon'].values, t['Lat'].values)      # (lat from -89.5, lon from -179.5)
        break
    if frac is None:
        print('LPJ-GUESS: no native cpool.out, 1 degree totals include whole coast cells (about 12 % high)')
    for d in sorted(glob.glob(f'{ROOT}/outdata/lpj_guess/*-*')):
        y0 = int(os.path.basename(d)[:4])
        try:
            cl = xr.open_dataset(f'{d}/cpool_remapcon_r360x180.nc', decode_times=False)
            cf = xr.open_dataset(f'{d}/cflux_remapcon_r360x180.nc', decode_times=False)
        except OSError:
            continue
        la = np.deg2rad(cl['lat'].values)
        ca = (RE ** 2 * np.deg2rad(1.0) * (np.sin(la + np.deg2rad(0.5)) - np.sin(la - np.deg2rad(0.5))))[:, None]
        if frac is not None:
            fr = frac[np.round(cl['lat'].values + 89.5).astype(int)][:, np.round(((cl['lon'].values + 180.0) % 360.0) - 0.5).astype(int)]
            ca = ca * fr
        for i in range(cl.sizes['time']):
            rows[y0 + i] = dict(
                land_c_PgC=float(np.nansum(cl['Total'].values[i] * ca)) / 1e12,
                veg_c_PgC=float(np.nansum(cl['VegC'].values[i] * ca)) / 1e12,
                soil_c_PgC=float(np.nansum((cl['SoilC'].values[i] + cl['LitterC'].values[i]) * ca)) / 1e12,
                nee_PgC=float(np.nansum(cf['NEE'].values[i] * ca)) / 1e12)
    return rows


if __name__ == '__main__':
    with Pool(NP) as p:
        rows = p.map(year, range(Y0, Y1 + 1))
    L = lpjg()
    for r in rows:
        r.update(L.get(r['year'], {}))
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False, float_format='%.7g')
    print(df.T.to_string())
