"""Ice-shelf basal melt by shelf from FESOM's fw under cavity nodes, against Rignot et al. 2013.

fw (water_flux, positive out of the ocean) under cavity nodes (ulevels_nod2D > 1) is minus the
basal melt.  Nodes are grouped into shelves by longitude/latitude boxes (coarse; small shelves
fall into 'other' sectors).  For each shelf: model melt in Gt/yr, its area-mean rate in m/yr,
the mean thermal forcing of the top wet layer (T minus the in-situ freezing point at that
depth, from temp/salt, annual) and Rignot 2013 (R13) and Adusumilli 2020 (A20, steady-state
2010-2018 where available) reference totals.  A shelf with modest thermal forcing but high
melt points at the transfer coefficient; high thermal forcing points at warm water on the
shelf (a circulation problem the coefficient cannot fix).

Usage:  ROOT=... ARM=PI200 Y0=1585 Y1=1599 MESH=<mesh diag> python3 scripts/analysis/cavity_melt_by_shelf.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ['ROOT']; ARM = os.environ['ARM']; MESH = os.environ['MESH']; Y0, Y1 = int(os.environ['Y0']), int(os.environ['Y1'])
# boxes: name, lon0, lon1 (degrees east, may wrap), lat0, lat1, R13 Gt/yr, A20 Gt/yr (approx.)
SHELVES = [('Filchner-Ronne', 275, 330, -84, -74.5, 155, 50), ('Ross', 158, 212, -86, -77, 48, 34), ('Amery', 62, 76, -74, -68, 36, 40),
           ('Larsen C/B', 293, 302, -70, -64, 21, 25), ('Fimbul-Riiser (0-30E, 330-360E)', 330, 30, -75, -68, 100, 80),
           ('Getz-Dotson-Crosson (230-250E)', 230, 250, -76, -73, 220, 250), ('Pine Island-Thwaites (250-262E)', 250, 262, -76, -74, 200, 190),
           ('Abbot-Venable-Bellingsh. (262-275E)', 262, 275, -74, -70, 80, 60), ('George VI-Wilkins-Stange (275-293E)', 275, 293, -74, -70, 150, 130),
           ('Totten-Moscow Univ. (112-125E)', 112, 125, -68, -66, 90, 80), ('Shackleton-West (80-112E)', 80, 112, -68, -65, 100, 90),
           ('Cook-Mertz-Ninnis (140-158E)', 140, 158, -70, -66, 30, 25), ('30-62E', 30, 62, -72, -66, 60, 50), ('212-230E (Sulzberger, Land)', 212, 230, -79, -73, 30, 25)]
with xr.open_dataset(MESH) as m:
    lat = m['lat'].values; lon = m['lon'].values % 360; ul = m['ulevels_nod2D'].values; z = m['nz1'].values
    area = m['nod_area'].values[ul - 1, np.arange(len(lat))]
cav = ul > 1


def ann(v):
    acc = 0; n = 0
    for y in range(Y0, Y1 + 1):
        f = f'{R}/{ARM}/outdata/fesom/{v}.fesom.{y}.nc'
        if not os.path.exists(f): continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = np.squeeze(d[v].values).astype('f8')
        a[np.abs(a) > 1e10] = np.nan; acc = acc + np.nanmean(a, 0); n += 1
    return acc / n


fw = ann('fw'); T = ann('temp'); S = ann('salt')
if T.shape[0] != len(lat): T = T.T; S = S.T
idx = np.arange(len(lat)); Tt = T[idx, ul - 1]; St = S[idx, ul - 1]; zt = z[ul - 1]
Tf = -0.0575 * St + 0.0901 - 7.61e-4 * (-abs(zt))          # Foldvik & Kvinge, as in cavity_param.F90 (in situ, approx.)
TF = Tt - Tf
GT = 1000.0 * 365.25 * 86400 / 1e12
melt = -fw * area * GT                                       # Gt/yr per node, positive = melting
print(f'{ARM} {Y0}-{Y1}: basal melt by shelf [Gt/yr], rate [m/yr], thermal forcing of the top wet layer [K]')
print(f'  {"shelf":<40}{"area Mkm2":>10}{"model":>8}{"R13":>6}{"A20":>6}{"rate":>7}{"TF":>7}')
used = np.zeros(len(lat), bool); tot = 0
for name, lo0, lo1, la0, la1, r13, a20 in SHELVES:
    inlon = (lon >= lo0) & (lon < lo1) if lo0 < lo1 else (lon >= lo0) | (lon < lo1)
    k = cav & inlon & (lat >= la0) & (lat < la1) & ~used; used |= k
    if not k.any(): print(f'  {name:<40}{"no cavity nodes":>20}'); continue
    mm = melt[k].sum(); tot += mm
    print(f'  {name:<40}{area[k].sum() / 1e12:10.3f}{mm:8.0f}{r13:6.0f}{a20:6.0f}{mm / GT / area[k].sum() * 365.25 * 86400 / 1000 * 1000 / 1000:7.2f}{np.average(TF[k], weights=area[k]):7.2f}')
k = cav & ~used
print(f'  {"other cavities":<40}{area[k].sum() / 1e12:10.3f}{melt[k].sum():8.0f}{"":>12}{np.average(TF[k], weights=area[k]) if k.any() else 0:14.2f}')
print(f'  {"TOTAL":<40}{area[cav].sum() / 1e12:10.3f}{melt[cav].sum():8.0f}{1325:6.0f}{1100:6.0f}{"":>7}{np.average(TF[cav], weights=area[cav]):7.2f}')
