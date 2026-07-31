"""Is the surface-albedo bias boreal-specific, or global across all land types?

The Siberian JJA surface SW deficit of 12.5 W/m2 splits almost evenly:
  ~7.0 W/m2  too little SW reaching the surface (excess cloud)
  ~5.5 W/m2  surface albedo too high (0.1753 model vs 0.1461 CERES)

The albedo half is independent of cloud and has never been touched. But before fixing it
for boreal needleleaf only, it must be established that the error IS boreal-specific. If
every land type is equally too bright, a vegetation-type-indexed fix would be papering over
a systematic error in the albedo scheme or in the prescribed albedo climatology, and would
mis-attribute a global problem to the boreal.

Albedo = SW_up / SW_down, model and CERES on the identical grid, box and land mask.
Regions are chosen to isolate distinct vegetation types and to AVOID permanent ice
(Greenland, Antarctica), whose albedo would swamp any box containing it. High-latitude
regions are evaluated in JJA so the comparison is snow-free.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
OBS = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
RUN = 'amip_pi_base'
Y0, Y1 = 1872, 1915
ACC = 3600.0
JJA, DJF, ANN = [5, 6, 7], [11, 0, 1], list(range(12))

# (name, lat range, lon range, months, expected dominant cover)
REGIONS = [
    ('Siberia boreal',     (55, 75), (60, 180),    JJA, 'needleleaf 3/4'),
    ('Canada boreal',      (55, 68), (-135, -70),  JJA, 'needleleaf 3/4'),
    ('Fennoscandia',       (58, 68), (10, 50),     JJA, 'needleleaf 3/4'),
    ('Siberian tundra',    (70, 77), (60, 160),    JJA, 'tundra 9'),
    ('Europe temperate',   (45, 55), (0, 30),      JJA, 'crops/broadleaf'),
    ('US Great Plains',    (35, 45), (-105, -95),  JJA, 'crops/grass'),
    ('Central Asia steppe',(42, 52), (55, 80),     JJA, 'short grass 2'),
    ('Sahara',             (18, 28), (-5, 30),     ANN, 'desert 8'),
    ('Amazon',             (-10, 0), (-70, -52),   ANN, 'evergreen broadleaf 6'),
    ('Congo',              (-5, 5),  (14, 28),     ANN, 'evergreen broadleaf 6'),
    ('India monsoon',      (12, 25), (74, 86),     JJA, 'crops 1'),
    ('Australia interior', (-30, -20), (120, 145), ANN, 'semidesert 11'),
    ('GLOBAL land 60S-70N',(-60, 70), (-180, 180), ANN, 'all'),
]

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values
_d = xr.open_dataset(f'{RT}/{RUN}/outdata/oifs/atm_remapped_1m_ssr_1m_{Y0}-{Y0}.nc')
lat, lon = _d['ssr'].lat.values, _d['ssr'].lon.values
_d.close()


def ceres(v):
    ds = xr.open_dataset(OBS)
    ds = xr.concat([ds.isel(lon=[-1]).assign_coords(lon=ds.lon.values[-1:] - 360.), ds,
                    ds.isel(lon=[0]).assign_coords(lon=ds.lon.values[:1] + 360.)], dim='lon')
    return ds[v].interp(lat=xr.DataArray(np.clip(lat, -89.5, 89.5), dims='y'),
                        lon=xr.DataArray(np.where(lon < 0, lon + 360, lon), dims='x')).values


def model(v):
    acc = []
    for y in range(Y0, Y1 + 1):
        f = f'{RT}/{RUN}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc'
        d = xr.open_dataset(f)
        acc.append(d[v].values / ACC)
        d.close()
    return np.mean(acc, axis=0)


def bm(a2d, la, lo):
    yi = (lat >= la[0]) & (lat <= la[1])
    l180 = ((lon + 180) % 360) - 180
    xi = (l180 >= lo[0]) & (l180 <= lo[1])
    s = a2d[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    L = lsm[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    m = L > 0.5
    if not m.any():
        return np.nan
    w = np.broadcast_to(np.cos(np.deg2rad(lat[yi]))[:, None], s.shape)
    return np.average(s[m], weights=w[m])


m_net, m_dn = model('ssr'), model('ssrd')
c_up, c_dn = ceres('sfc_sw_up_all_clim'), ceres('sfc_sw_down_all_clim')

print(f'Surface albedo, model ({RUN}, {Y0}-{Y1}) vs CERES, identical box + land mask\n')
print(f"  {'region':22s}{'cover':18s}{'model':>8s}{'CERES':>8s}{'diff':>8s}{'W/m2':>8s}")
for name, la, lo, mo, cover in REGIONS:
    mn = np.mean([bm(m_net[k], la, lo) for k in mo])
    md = np.mean([bm(m_dn[k], la, lo) for k in mo])
    cu = np.mean([bm(c_up[k], la, lo) for k in mo])
    cd = np.mean([bm(c_dn[k], la, lo) for k in mo])
    if not np.isfinite(mn) or not np.isfinite(cu):
        print(f'  {name:22s}{cover:18s}      no land points')
        continue
    am, ao = (md - mn) / md, cu / cd
    # W/m2 of absorbed SW lost to the albedo error alone, at the MODEL's own SW_down
    print(f'  {name:22s}{cover:18s}{am:8.4f}{ao:8.4f}{am-ao:+8.4f}{-(am-ao)*md:+8.2f}')
print('\n  diff > 0 means the model is TOO BRIGHT; W/m2 is the absorbed SW lost to that,')
print('  evaluated at the model\'s own SW_down so it is independent of the cloud error.')
