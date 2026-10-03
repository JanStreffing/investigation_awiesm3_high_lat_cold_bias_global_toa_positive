"""Southern Ocean JJA warm bias over sea ice: missing ice, or ice that is too warm?

The JJA T2m bias against ERA5 peaks at 60-75S right over the winter pack.  Two readings
need different levers:
  EXTENT   the model has open water (about -1.8 C) where the observed pack sits at -15 C,
           so the bias is the ice edge, and the lever is ocean heat / ice growth;
  PACK     the model has ice, but the air above it is too warm, so the lever is the
           surface energy balance over ice: downward LW from supercooled cloud, the
           turbulent coupling of a stable boundary layer, or conduction through thin ice.

Classes on the JJA-mean concentration, model ci (OpenIFS remapped) against HadISST sic,
55-80S ocean, on the 1x1 CERES/HadISST grid (model nearest neighbour):
  pack     both >= 0.8           loose    obs >= 0.8, model 0.15-0.8
  missing  obs >= 0.15, model < 0.15      extra    model >= 0.15, obs < 0.15
  open     both < 0.15
For each class: area, T2m and skin bias against ERA5 (1990-2014), and its contribution to
the 55-80S ocean mean bias.  In the pack class also the surface energy terms, model
against CERES EBAF surface fluxes (a radiative-transfer product, weakest over ice in polar
night, so a guide rather than a verdict).  IFS fluxes are accumulated J/m2 per output hour
(/3600), positive downward.  Two windows are shown so a change in the pack terms between
them can be read as the effect of the cloud-scheme change plus the coupled cooling.

Usage:  ARM=PI200 WINDOWS=1480-1499,1520-1539 python3 scripts/analysis/so_winter_ice_warm_bias.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); ACCF = float(os.environ.get('ACC', 3600.0))   # ACC=1 for runs whose fluxes are already W/m2 (awicm3 v3.3)
WINDOWS = [tuple(map(int, w.split('-'))) for w in os.environ.get('WINDOWS', '1480-1499,1520-1539').split(',')]
CER = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
SIC = '/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc'
E5 = '/work/ab0246/a270092/obs/era5/netcdf'
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'
JJA = [5, 6, 7]

with xr.open_dataset(CER, decode_times=False) as dc:
    clat = np.squeeze(dc['lat'].values); clon = np.squeeze(dc['lon'].values)
    C = {k: np.squeeze(dc[v].values)[JJA].mean(0) for k, v in (
        ('LWdn', 'sfc_lw_down_all_clim'), ('LWdn_clr', 'sfc_lw_down_clr_c_clim'),
        ('LWup', 'sfc_lw_up_all_clim'), ('SWdn', 'sfc_sw_down_all_clim'), ('SWup', 'sfc_sw_up_all_clim'))}
with xr.open_dataset(SIC, decode_times=False) as ds:
    sic = np.squeeze(ds['sic'].values); slat = np.squeeze(ds['latitude'].values)
if np.nanmax(sic) > 1.5: sic = sic / 100.0
if slat[0] > slat[-1]: sic = sic[:, ::-1, :]
sic = np.nan_to_num(sic[JJA].mean(0), nan=0.0)                 # HadISST lon is already 0.5..359.5


def era5(f, v):
    with xr.open_dataset(f'{E5}/{f}') as d:
        c = d[v].sel(time=d['time.year'].isin(range(1990, 2015)))
        c = c.sel(time=c['time.month'].isin([6, 7, 8])).mean('time')
        c = c.assign_coords(longitude=c['longitude'] % 360).sortby('longitude').sortby('latitude')
        return np.asarray(c.interp(latitude=('y', clat), longitude=('x', clon)).values, float)


O = {'T2m': era5('T2M.nc', 't2m'), 'skt': era5('SKT_mon.nc', 'skt')}
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = lsm[0] if lsm.ndim == 3 else lsm
    mlat = np.squeeze(d['lat'].values); mlon = np.squeeze(d['lon'].values)
iy = np.abs(clat[:, None] - mlat[None, :]).argmin(1)
ix = np.abs(((clon[:, None] - mlon[None, :] + 180) % 360 - 180)).argmin(1)
rg = lambda a: a[..., iy, :][..., ix]
LAT = np.broadcast_to(clat[:, None], sic.shape)
W = np.cos(np.deg2rad(LAT))
OCN = (rg(lsm) <= 0.5) & (LAT >= -80) & (LAT <= -55)


def mload(v, y):
    for p in (f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc',
              f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc'):
        if os.path.exists(p):
            with xr.open_dataset(p, decode_times=False) as d:
                k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
                return np.squeeze(d[k].values).astype('f8')
    if v == 'strdc':                                   # older runs (awicm3 v3.3) have no clear-sky LW down
        return None
    raise FileNotFoundError(f'{v} {y}')


VARS = ('2t', 'skt', 'ci', 'strd', 'strdc', 'str', 'ssr', 'sshf', 'slhf', 'tclw', 'tciw', 'tcc')
ACC = {'strd', 'strdc', 'str', 'ssr', 'sshf', 'slhf'}
for y0, y1 in WINDOWS:
    M = {v: 0.0 for v in VARS}
    for y in range(y0, y1 + 1):
        for v in VARS:
            a = mload(v, y)
            if a is None:
                M[v] = None; continue
            a = a[JJA].mean(0)
            M[v] = M[v] + rg(a / ACCF if v in ACC else a) / (y1 - y0 + 1)
    ci = np.nan_to_num(M['ci'], nan=0.0)
    K = {'pack': (ci >= 0.8) & (sic >= 0.8), 'loose': (sic >= 0.8) & (ci >= 0.15) & (ci < 0.8),
         'missing': (sic >= 0.15) & (ci < 0.15), 'extra': (ci >= 0.15) & (sic < 0.15),
         'open': (ci < 0.15) & (sic < 0.15)}
    K['obs 0.15-0.8, model>=0.15'] = (sic >= 0.15) & (sic < 0.8) & (ci >= 0.15)
    bT = M['2t'] - O['T2m']; bS = M['skt'] - O['skt']
    ok = OCN & np.isfinite(bT)
    wt = W[ok].sum(); area = lambda k: float((W * 111.195e3 ** 2)[k].sum() / 1e12)
    print(f'\n{ARM} {y0}-{y1}  JJA, 55-80S ocean: T2m minus ERA5 {np.average(bT[ok], weights=W[ok]):+.2f} K'
          f'   model ice area-weighted extent {area(ok & (ci >= 0.15)):.2f}, HadISST {area(ok & (sic >= 0.15)):.2f} M km2 (1x1 cells, 55-80S)')
    print(f'  {"class":<28}{"M km2":>7}{"T2m bias":>10}{"skin bias":>10}{"contrib":>9}{"model ci":>9}{"obs sic":>8}')
    for n, k in K.items():
        k = k & ok
        if not k.any(): continue
        av = lambda f: np.average(f[k], weights=W[k])
        print(f'  {n:<28}{area(k):7.2f}{av(bT):+10.2f}{av(bS):+10.2f}{(bT * W)[k].sum() / wt:+9.2f}{av(ci):9.2f}{av(sic):8.2f}')
    k = K['pack'] & ok
    av = lambda f: np.average(f[k], weights=W[k])
    print('  pack, surface terms [W/m2, positive down]:')
    print(f'    LW down all   model {av(M["strd"]):7.1f}   CERES {av(C["LWdn"]):7.1f}')
    if M['strdc'] is not None:
        print(f'    LW down clear model {av(M["strdc"]):7.1f}   CERES {av(C["LWdn_clr"]):7.1f}')
        print(f'    LW cloud eff. model {av(M["strd"] - M["strdc"]):7.1f}   CERES {av(C["LWdn"] - C["LWdn_clr"]):7.1f}')
    print(f'    LW net        model {av(M["str"]):7.1f}   CERES {av(C["LWdn"] - C["LWup"]):7.1f}')
    print(f'    SW net        model {av(M["ssr"]):7.1f}   CERES {av(C["SWdn"] - C["SWup"]):7.1f}')
    print(f'    sensible {av(M["sshf"]):6.1f}   latent {av(M["slhf"]):6.1f}   net {av(M["str"] + M["ssr"] + M["sshf"] + M["slhf"]):6.1f}')
    print(f'    T2m - skin model {av(M["2t"] - M["skt"]):+.2f} K, ERA5 {av(O["T2m"] - O["skt"]):+.2f} K'
          f'   tcc {av(M["tcc"]):.2f}  tclw {av(M["tclw"]) * 1e3:.1f} g/m2  tciw {av(M["tciw"]) * 1e3:.1f} g/m2')
