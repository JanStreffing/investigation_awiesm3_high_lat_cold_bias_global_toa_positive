"""v3.5.0-cc spin-up (1875-1884) against the tuning-line iceberg run PICAL_crunveg_icb (2195-2209).

By latitude band: T2m (all, land, ocean), net TOA solar, forecast albedo and snow depth over land;
September area of MLD2 > 1000 m south of 50S; LPJ-GUESS carbon pools per band (kgC/m2, land cells) at the
start and end of the spin-up and in the tuning line; the land CO2 flux as OpenIFS has it (co2nbf) against
LPJ-GUESS NEE.

Usage: spinup_v350_vs_tuning_line.py   (module analysis-toolbox/python-04.2026)
"""
import os, glob
import numpy as np, pandas as pd, xarray as xr
S = '/albedo/work/projects/p_awiesm3_cmip7/lahajd001/runtime/awiesm3-v3.5.0-cc/AWI-ESM3-CC-LR-CMIP7-Spinup'
R = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4/PICAL_crunveg_icb'
YS, YR_ = range(1875, 1885), range(2195, 2210)
BANDS = [(-90, -60), (-60, -30), (-30, 30), (30, 60), (60, 90)]
so = lambda v, y: xr.open_dataset(f'{S}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc')[v]
ro = lambda v, y: xr.open_dataset(f'{R}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc')[v]
def clim(f, v, ys):
    return sum(f(v, y).mean('time_counter').load() for y in ys) / len(ys)
lsm = so('lsm', 1880).isel(time_counter=0).load()
w = np.cos(np.deg2rad(lsm['lat'])) + 0 * lsm
def band(x, a, b, m=None):
    k = (x['lat'] >= a) & (x['lat'] < b)
    ww = w.where(k, 0) if m is None else w.where(k & m, 0)
    return float((x * ww).sum() / ww.sum())
print('band          T2m s    ref   diff | land diff  ocean diff | TOA solar diff W/m2 | land albedo s  ref | land snow depth s  ref (m w.e.)')
F = {v: (clim(so, v, YS), clim(ro, v, YR_)) for v in ('2t', 'tsr', 'fal', 'sd')}
for a, b in BANDS + [(-90, 90)]:
    t = [band(F['2t'][i], a, b) - 273.15 for i in (0, 1)]
    dl = band(F['2t'][0] - F['2t'][1], a, b, lsm > 0.5); do = band(F['2t'][0] - F['2t'][1], a, b, lsm <= 0.5)
    ds = band(F['tsr'][0] - F['tsr'][1], a, b) / 3600
    al = [band(F['fal'][i], a, b, lsm > 0.5) for i in (0, 1)]; sd = [band(F['sd'][i], a, b, lsm > 0.5) for i in (0, 1)]
    print(f'{a:4d}..{b:3d}  {t[0]:7.2f} {t[1]:7.2f} {t[0]-t[1]:+6.2f} |  {dl:+6.2f}    {do:+6.2f}   |      {ds:+6.2f}         |   {al[0]:.3f} {al[1]:.3f}   |   {sd[0]:.3f} {sd[1]:.3f}')
# seasonal T2m difference 60-90N and 30-60N land
d = sum(so('2t', y).load().values for y in YS) / len(YS) - sum(ro('2t', y).load().values for y in YR_) / len(YR_)
for a, b in ((60, 90), (30, 60)):
    k = ((lsm['lat'] >= a) & (lsm['lat'] < b)).values[:, None] & (lsm.values > 0.5); ww = w.values * k
    print(f'land {a}-{b}N T2m diff by month:', ' '.join(f'{(d[m]*ww).sum()/ww.sum():+.1f}' for m in range(12)))
# Southern Ocean deep mixed layers
md = xr.open_dataset('/albedo/work/projects/p_awiesm3_cmip7/jstreffi/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lat = md['lat'].values; area = md['nod_area'].values[0]; oc = md['ulevels_nod2D'].values <= 1
def mld(root, ys):
    return np.mean([area[(lat < -50) & oc & (np.abs(xr.open_dataset(f'{root}/outdata/fesom/MLD2.fesom.{y}.nc')['MLD2'].isel(time=8).values) > 1000)].sum() / 1e12 for y in ys])
print(f'SO September MLD2>1000 m area: spin-up {mld(S, YS):.2f}, tuning line {mld(R, YR_):.2f} Mkm2')
# LPJ-GUESS pools per band
def sp(d, i):
    c = xr.open_dataset(f'{S}/outdata/lpj_guess/{d}/cpool_remapcon_r360x180.nc', decode_times=False)
    la = c['lat'].values; ca = np.cos(np.deg2rad(la))[:, None] + 0 * c['VegC'].values[i]
    out = []
    for a, b in BANDS[1:]:
        k = (la >= a) & (la < b); r = []
        for v in ('VegC', 'SoilC', 'LitterC'):
            x = c[v].values[i][k]; ww = np.where(np.isfinite(x), ca[k], 0)
            r.append(np.nansum(x * ww) / ww.sum())
        out.append(r)
    return out
def rf(d, yr):
    t = pd.read_csv(f'{R}/outdata/lpj_guess/{d}/cpool.out', sep=r'\s+'); t = t[t.Year == yr]
    out = []
    for a, b in BANDS[1:]:
        q = t[(t.Lat >= a) & (t.Lat < b)]
        out.append([q.VegC.mean(), q.SoilC.mean(), q.LitterC.mean()])
    return out
A, B_, C = sp('18500101-18541231', 0), sp('18800101-18841231', 4), rf('22000101-22091231', 2205)
print('LPJ-GUESS kgC/m2 over land cells (VegC SoilC LitterC):  spin-up 1850 | spin-up 1884 | tuning line 2205')
for (a, b), x, y, z in zip(BANDS[1:], A, B_, C):
    print(f'{a:4d}..{b:3d}   ' + ' | '.join(' '.join(f'{v:6.2f}' for v in r) for r in (x, y, z)))
tf = lambda d, i: xr.open_dataset(f'{S}/outdata/lpj_guess/{d}/treeFrac_yearly_remapcon_r360x180.nc', decode_times=False)
for d, i, lab in (('18500101-18541231', 0, '1850'), ('18800101-18841231', 4, '1884')):
    t = tf(d, i); v = list(t.data_vars)[0]; x = t[v].values[i]; la = t['lat'].values; ca = np.cos(np.deg2rad(la))[:, None] + 0 * x
    print(f'tree fraction {lab} ({t[v].attrs.get("units")}):', ' '.join(f'{a}..{b}: {np.nansum((x*ca)[(la>=a)&(la<b)])/np.where(np.isfinite(x),ca,0)[(la>=a)&(la<b)].sum():.3f}' for a, b in BANDS[1:]))
# co2nbf as OpenIFS has it
x = sum(so('co2nbf', y).load() for y in YS) / len(YS)
AG = 4 * np.pi * 6.371e6 ** 2; K = AG * 365.25 * 86400 * 12 / 44 / 1e12
wl = w.where(lsm > 0.5, 0); wo = w.where(lsm <= 0.5, 0)
print('co2nbf PgC/yr: land cells', float((x.mean('time_counter') * wl).sum() / w.sum()) * K, ' ocean cells', float((x.mean('time_counter') * wo).sum() / w.sum()) * K)
print('co2nbf land, by month (PgC/yr):', ' '.join(f'{float((x.isel(time_counter=m) * wl).sum() / w.sum()) * K:.0f}' for m in range(12)))
nh = wl.where(lsm['lat'] > 30, 0)
print('co2nbf land north of 30N, by month:', ' '.join(f'{float((x.isel(time_counter=m) * nh).sum() / w.sum()) * K:.0f}' for m in range(12)))
print('co2nbf min/max kg m-2 s-1:', float(x.min()), float(x.max()))
