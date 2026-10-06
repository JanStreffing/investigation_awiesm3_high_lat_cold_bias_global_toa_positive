"""Southern Ocean westerlies and wind stress: AMIP and coupled runs against ERA5, like for like.

Zonal mean over the ocean of the 10 m zonal wind (10u) and of the eastward turbulent surface stress
(ewss, accumulated N m-2 s per hour, divided by 3600), annual mean, for runs on this machine, against
ERA5 1990-2014 (u10 from obs/era5/netcdf/uas_vas.nc; stress from the DKRZ pool, E5 sf fc 1M param 180).
ERA5 is brought to the model grid (bilinear) and masked with the model's land-sea mask.
Reports jet maximum, its latitude and the 50-65S mean for both quantities, and the ratio
stress / u10^2 over 50-65S as a bulk drag measure.

Usage (levante): python so_wind_stress_amip.py
"""
import os, glob, subprocess, sys
import numpy as np
import xarray as xr

OUT = os.path.dirname(os.path.abspath(__file__))
RUNS = [('AMIP present day, v3.5 (amip_PD1_v35)', '/work/bb1469/a270092/runtime/oifsamip-cy48/amip_PD1_v35/outdata/oifs', 1990, 2014),
        ('AMIP present day, control (amip_PD0_ctl)', '/work/bb1469/a270092/runtime/oifsamip-cy48/amip_PD0_ctl/outdata/oifs', 1990, 2014),
        ('AMIP pre-industrial SST (amip_picontrol)', '/work/bb1469/a270092/runtime/oifsamip-cy48/amip_picontrol/outdata/oifs', 1880, 1917),
        ('coupled PI, PICAL_ccnice 2085-2099', '/work/bb1469/a270092/runtime/awiesm3-v3.4/PICAL_ccnice/outdata/oifs', 2085, 2099)]
E5U = '/work/ab0246/a270092/obs/era5/netcdf/uas_vas.nc'
E5S = '/pool/data/ERA5/E5/sf/fc/1M/180'


def files(d, var, y):
    f = glob.glob(f'{d}/atm_remapped_1m_{var}_{y}-{y}.nc') + glob.glob(f'{d}/atm_remapped_1m_{var}_1m_{y}-{y}.nc')
    return f[0] if f else None


def tmean(d, var, y0, y1):
    acc, n = 0, 0
    for y in range(y0, y1 + 1):
        f = files(d, var, y)
        if f is None:
            continue
        with xr.open_dataset(f, decode_times=False) as ds:
            k = [c for c in ds.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            x = ds[k].mean([dd for dd in ds[k].dims if dd not in ('lat', 'lon')])
            acc = acc + x.load(); n += 1
    return (acc / n).sortby('lat'), n


d0 = RUNS[0][1]
grid = files(d0, '10u', 1990)
with xr.open_dataset(files(d0, 'lsm', 1990) or '/work/bb1469/a270092/runtime/awiesm3-v3.4/PICAL_ccnice/outdata/oifs/atm_remapped_1m_lsm_2090-2090.nc', decode_times=False) as ds:
    if 'lsm' in ds:
        lsm = ds['lsm'].mean([dd for dd in ds['lsm'].dims if dd not in ('lat', 'lon')]).sortby('lat').load()
    else:
        lsm = None
if lsm is None:
    sys.exit('no lsm file found next to the AMIP output')
ocean = lsm < 0.5

# ERA5 stress on the model grid
es = f'{OUT}/era5_ewss_1990-2014_modelgrid.nc'
if not os.path.exists(es):
    fl = ' '.join(f'{E5S}/E5sf12_1M_{y}_180.grb' for y in range(1990, 2015))
    subprocess.run(f'cdo -s -f nc -b F32 -remapbil,{grid} -setgridtype,regular -timmean -mergetime {fl} {es}', shell=True, check=True)
with xr.open_dataset(es, decode_times=False) as ds:
    k = [c for c in ds.data_vars if 'lat' in ds[c].dims][0]
    e_s = ds[k].squeeze().sortby('lat').load()
with xr.open_dataset(E5U) as ds:
    e_u = ds['u10'].sel(time=ds['time.year'].isin(range(1990, 2015))).mean('time').rename({'latitude': 'lat', 'longitude': 'lon'}).sortby('lat').load()
e_u = e_u.interp(lat=lsm['lat'], lon=lsm['lon'])
# ERA5 monthly means of accumulations: find the divisor that gives a plausible stress (N/m2)
es_raw = float(np.abs(e_s.where(ocean).mean('lon')).max())
div = 1.0 if es_raw < 1 else (3600.0 if es_raw < 3600 else 86400.0)
e_s = e_s / div


def stats(u, s):
    zu = u.where(ocean).mean('lon').sel(lat=slice(-75, -30)); zs = s.where(ocean).mean('lon').sel(lat=slice(-75, -30))
    w = np.cos(np.deg2rad(zu['lat'])); b = (zu['lat'] >= -65) & (zu['lat'] <= -50)
    bm = lambda z: float((z * w).where(b).sum() / w.where(b).sum())
    return (float(zu.max()), float(zu['lat'][int(zu.argmax())]), bm(zu), float(zs.max()), float(zs['lat'][int(zs.argmax())]), bm(zs), bm(zs) / bm(zu) ** 2)


E = stats(e_u, e_s)
out = [f'Zonal mean over the ocean, annual. ERA5 1990-2014 (stress file divided by {div:g}).',
       f"{'':46s}{'u10 max':>8s}{'lat':>7s}{'u10 50-65S':>11s}{'tau max':>9s}{'lat':>7s}{'tau 50-65S':>11s}{'tau/u10^2':>11s}{'years':>6s}",
       f"{'ERA5':46s}{E[0]:8.2f}{E[1]:7.1f}{E[2]:11.2f}{E[3]:9.3f}{E[4]:7.1f}{E[5]:11.3f}{1e3 * E[6]:11.2f}"]
for name, d, y0, y1 in RUNS:
    u, n = tmean(d, '10u', y0, y1); s, _ = tmean(d, 'ewss', y0, y1)
    u = u.interp(lat=lsm['lat'], lon=lsm['lon']); s = s.interp(lat=lsm['lat'], lon=lsm['lon']) / 3600.0
    r = stats(u, s)
    out.append(f"{name:46s}{r[0]:8.2f}{r[1]:7.1f}{r[2]:11.2f}{r[3]:9.3f}{r[4]:7.1f}{r[5]:11.3f}{1e3 * r[6]:11.2f}{n:6d}")
    out.append(f"{'   ratio to ERA5':46s}{r[0] / E[0]:8.2f}{'':7s}{r[2] / E[2]:11.2f}{r[3] / E[3]:9.2f}{'':7s}{r[5] / E[5]:11.2f}{r[6] / E[6]:11.2f}")
txt = '\n'.join(out)
print(txt)
open(f'{OUT}/so_wind_stress_amip.txt', 'w').write(txt + '\n')
