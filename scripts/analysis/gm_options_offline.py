"""Offline GM coefficient for FESOM's scaling options, from a run's own monthly N2 and MLD1,
to see which physically based option keeps K high where Southern Ocean mode/intermediate
water forms and low where the Labrador Sea convects. Follows oce_fer_gm.F90:init_Redi_GM.

  horizontal  res:    min(1, (2 A / (100 km)^2)^0.5)            (scaling_resolution, order 2)
              rossby: Fermi(min(reso/Rd, 5); x0 1.5, sigma 0.15)  (scaling_rossby)
                      Rd = min(c1/max(|f|,1e-6), 200 km), c1 = max(0.1, sum(N dz)/pi), reso = 2 sqrt(A/pi)
  vertical    zexp:   0.6 + 0.4 exp(-|z|/500 m)                     (scaling_gmzexp, production)
              ferr:   clip(N2(z) / N2(first interface below MLD1), 0.2, 1)  (scaling_ferreira, K_GM_bvref 1)
Options reported (K_GM_max 2500 m2/s): res*zexp (production), res*ferr, rossby*zexp, rossby*ferr.
Monthly fields averaged over the years given; regions by node lon/lat; depth bands by interface depth.
Slope tapering (odm95) and the K_GM_min floor are left out, since they act on all options alike.
Usage: python gm_options_offline.py <mesh.diag.nc> <outdata_fesom> <y0> <y1>
"""
import sys, numpy as np
from netCDF4 import Dataset
diag, d, y0, y1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
KMAX = 2500.0
with Dataset(diag) as nc:
    A = np.asarray(nc.variables['nod_area'][0], 'f8')
    zb = np.abs(np.asarray(nc.variables['nz'][:], 'f8'))              # 48 interface depths
    lon = np.asarray(nc.variables['lon'][:], 'f8'); lat = np.asarray(nc.variables['lat'][:], 'f8')
if np.abs(lon).max() < 7: lon, lat = np.rad2deg(lon), np.rad2deg(lat)
lon = ((lon + 180) % 360) - 180
f = np.abs(2 * 7.292e-5 * np.sin(np.deg2rad(lat)))
reso = 2 * np.sqrt(A / np.pi)
s_res = np.minimum(1.0, np.sqrt(2 * A / 1e10))
dz = np.diff(zb)                                                       # 47 layer thicknesses
zexp = 0.6 + 0.4 * np.exp(-zb / 500.0)                                 # at interfaces
REG = {'Labrador 56-62N 60-50W': (lat >= 56) & (lat < 62) & (lon >= -60) & (lon < -50),
       'subpolar N Atl 45-66N':  (lat >= 45) & (lat < 66) & (lon >= -60) & (lon <= -10),
       'Nordic 66-80N 20W-20E':  (lat >= 66) & (lat < 80) & (lon >= -20) & (lon <= 20),
       'N Atl subtrop 20-40N':   (lat >= 20) & (lat < 40) & (lon >= -70) & (lon <= -20),
       'SO 40-55S (SAMW)':       (lat < -40) & (lat >= -55),
       'SO 55-65S':              (lat < -55) & (lat >= -65)}
BAND = {'0-200': zb < 200, '200-1000': (zb >= 200) & (zb < 1000), '1000-2000': (zb >= 1000) & (zb < 2000)}
acc = {m: [np.zeros((len(lat), 48)), np.zeros(len(lat)), 0] for m in (2, 8)}  # March, September
for y in range(y0, y1 + 1):
    with Dataset(f'{d}/N2.fesom.{y}.nc') as nc, Dataset(f'{d}/MLD1.fesom.{y}.nc') as nm:
        for m in acc:
            n2 = np.asarray(nc.variables['N2'][m], 'f8'); n2[~np.isfinite(n2) | (np.abs(n2) > 1)] = np.nan
            mld = np.abs(np.asarray(nm.variables['MLD1'][m], 'f8'))
            acc[m][0] += np.nan_to_num(n2); acc[m][1] += np.nan_to_num(mld); acc[m][2] += 1
for m, (n2s, mlds, n) in acc.items():
    n2 = n2s / n; mld = mlds / n
    wet = np.abs(n2) > 0                                               # interfaces with data
    N = np.sqrt(np.maximum(n2, 0.0))
    c1 = np.maximum(0.1, (0.5 * (N[:, :-1] + N[:, 1:]) * dz).sum(1) / np.pi)
    Rd = np.minimum(c1 / np.maximum(f, 1e-6), 2e5)
    s_ros = 1.0 / (1.0 + np.exp(-(np.minimum(reso / Rd, 5.0) - 1.5) / 0.15))
    ib = np.clip(np.searchsorted(zb, mld) , 1, 47)                     # first interface below MLD1
    bvref = np.maximum(n2[np.arange(len(lat)), ib], 1e-6)
    ferr = np.clip(n2 / bvref[:, None], 0.2, 1.0)
    OPT = {'res*zexp (prod)': s_res[:, None] * zexp[None, :], 'res*ferreira': s_res[:, None] * ferr,
           'rossby*zexp': s_ros[:, None] * zexp[None, :], 'rossby*ferreira': s_ros[:, None] * ferr}
    print(f'\n=== {"March" if m == 2 else "September"} (mean over {n} yrs), K in m2/s, area-weighted; MLD1 and Rd shown for context')
    print(f'{"region":24s} {"MLD1":>6s} {"Rd km":>6s} ' + ' '.join(f'{o:>18s}' for o in OPT))
    for r, rm in REG.items():
        row = f'{r:24s} {np.average(mld[rm], weights=A[rm]):6.0f} {np.average(Rd[rm], weights=A[rm])/1e3:6.1f} '
        cells = []
        for o, sc in OPT.items():
            vals = []
            for b, bm in BAND.items():
                k = KMAX * sc[np.ix_(rm, bm)]; w = (A[rm][:, None] * wet[np.ix_(rm, bm)])
                vals.append((k * w).sum() / max(w.sum(), 1e-30))
            cells.append('/'.join(f'{v:4.0f}' for v in vals))
        print(row + ' '.join(f'{c:>18s}' for c in cells))
print('\ncells: mean K over 0-200 / 200-1000 / 1000-2000 m')
