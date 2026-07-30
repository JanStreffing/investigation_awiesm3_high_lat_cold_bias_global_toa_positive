"""The period-clean measurement: same model, run in the reference period itself.

Every target in this campaign has been scored across a ~130-year gap -- runs at 1870s SST
with 1850 GHG, references being ERA5 1990-2014 and the CERES EBAF 07/2005-06/2015
climatology. Yesterday that gap was estimated INDIRECTLY: ERA5 back-periods chained through
HadCRUT5's global series and scaled by an amplification factor for T2m, and ERA5's
model-derived TOA fluxes for the energy target. Both were flagged as bounds, not
measurements.

amip_presentday removes the need for any of that. It is the same model and configuration
run over 1989-2015 with transient GHG, so:

  control (1872-1915)     vs ERA5 1990-2014  =  bias + period offset   <- what we tuned against
  presentday (1990-2014)  vs ERA5 1990-2014  =  bias alone             <- the honest number
  difference of the two                      =  the period offset, measured BY THE MODEL

That last line is the point. It replaces the HadCRUT5 chain with an internal, self-consistent
estimate on the identical grid, mask, box and season.

Same treatment for radiation against CERES: SO 45-65S SW CRE and global net TOA.

CAVEATS.
  * presentday uses transient historical GHG/aerosol while the control is fixed 1850, which
    is the intended contrast, but it also means aerosol forcing differs -- part of the
    offset is aerosol, not just greenhouse warming. This measures the TOTAL epoch
    difference, which is what we want for target-setting, not the CO2 contribution alone.
  * CERES climatology is 07/2005-06/2015 while presentday is 1990-2014, so a small residual
    period gap remains on the radiation side (~12 yr, not ~130).
  * ERA5 T2m is assimilated and trustworthy; ERA5 is not used here for radiation.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

ACC = 3600.0
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
OBS = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
ERA5 = '/work/ab0246/a270092/obs/era5/netcdf/T2M.nc'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
JJA = [5, 6, 7]
SIB = ((55, 75), (60, 180))
ARMS = [('control  (1872-1915, 1850 GHG)', 'amip_pi_base', 1872, 1915),
        ('presentday (1990-2014, transient)', 'amip_presentday', 1990, 2014)]

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def load(run, var, y0, y1):
    acc = []
    for y in range(y0, y1 + 1):
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            return None, None, None
        ds = xr.open_dataset(f)
        a = ds[var].values
        if var in {'tsr', 'tsrc', 'ttr', 'ssr'}:
            a = a / ACC
        lat, lon = ds[var].lat.values, ds[var].lon.values
        ds.close()
        acc.append(a)
    return np.mean(acc, axis=0), lat, lon


def sel(f, lat, lon, la, lo=(-180, 180), sfc='all', months=None):
    yi = (lat >= la[0]) & (lat <= la[1])
    l180 = ((lon + 180) % 360) - 180
    xi = (l180 >= lo[0]) & (l180 <= lo[1])
    mo = months if months is not None else list(range(12))
    sub = f[np.ix_(mo, np.where(yi)[0], np.where(xi)[0])]
    L = lsm[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    m = np.ones(L.shape, bool) if sfc == 'all' else (L > 0.5 if sfc == 'land' else L <= 0.5)
    w = np.broadcast_to(np.cos(np.deg2rad(lat[yi]))[:, None], sub.shape[1:])
    return float(np.mean([np.average(s[m], weights=w[m]) for s in sub]))


def era5_jja_siberia(lat, lon):
    """ERA5 1990-2014 JJA Siberian land mean, on the model grid with the model mask."""
    ds = xr.open_dataset(ERA5)
    da = ds['t2m']
    la = 'latitude' if 'latitude' in da.dims else 'lat'
    lo = 'longitude' if 'longitude' in da.dims else 'lon'
    da = da.rename({la: 'lat', lo: 'lon'}).sortby('lat')
    da = da.groupby('time.month').mean('time')
    da = xr.concat([da.isel(lon=[-1]).assign_coords(lon=da.lon.values[-1:] - 360.0), da,
                    da.isel(lon=[0]).assign_coords(lon=da.lon.values[:1] + 360.0)], dim='lon')
    tl, tlo = np.clip(lat, -89.5, 89.5), np.where(lon < 0, lon + 360, lon)
    g = da.interp(lat=xr.DataArray(tl, dims='y'), lon=xr.DataArray(tlo, dims='x')).values
    return sel(g, lat, lon, SIB[0], SIB[1], 'land', JJA) - 273.15


def ceres(lat, lon, v):
    ds = xr.open_dataset(OBS)
    ds = xr.concat([ds.isel(lon=[-1]).assign_coords(lon=ds.lon.values[-1:] - 360.0), ds,
                    ds.isel(lon=[0]).assign_coords(lon=ds.lon.values[:1] + 360.0)], dim='lon')
    tl, tlo = np.clip(lat, -89.5, 89.5), np.where(lon < 0, lon + 360, lon)
    return ds[v].interp(lat=xr.DataArray(tl, dims='y'), lon=xr.DataArray(tlo, dims='x')).values


R = {}
for name, run, y0, y1 in ARMS:
    t2m, lat, lon = load(run, '2t', y0, y1)
    if t2m is None:
        print(f'  !! {name}: incomplete'); continue
    tsr, _, _ = load(run, 'tsr', y0, y1)
    tsrc, _, _ = load(run, 'tsrc', y0, y1)
    ttr, _, _ = load(run, 'ttr', y0, y1)
    R[name] = dict(
        sib=sel(t2m, lat, lon, SIB[0], SIB[1], 'land', JJA) - 273.15,
        so_cre=sel(tsr - tsrc, lat, lon, (-65, -45), sfc='ocean'),
        toa=sel(tsr + ttr, lat, lon, (-90, 90)))
    R[name]['_grid'] = (lat, lon)

lat, lon = list(R.values())[0]['_grid']
E5 = era5_jja_siberia(lat, lon)
C_cre = sel(ceres(lat, lon, 'toa_cre_sw_clim'), lat, lon, (-65, -45), sfc='ocean')
C_toa = sel(ceres(lat, lon, 'toa_net_all_clim'), lat, lon, (-90, 90))

print(f'Observational references: ERA5 1990-2014 Siberian JJA T2m = {E5:.2f} C')
print(f'                          CERES 2005-2015 SO SW CRE = {C_cre:.2f}, '
      f'global net TOA = {C_toa:.2f} W/m2\n')
print(f"  {'arm':36s} {'Sib JJA T2m':>12s} {'bias':>8s} | {'SO SW CRE':>10s} {'bias':>8s} | "
      f"{'net TOA':>8s} {'bias':>8s}")
for k, v in R.items():
    print(f"  {k:36s} {v['sib']:12.2f} {v['sib']-E5:+8.2f} | {v['so_cre']:10.2f} "
          f"{v['so_cre']-C_cre:+8.2f} | {v['toa']:8.2f} {v['toa']-C_toa:+8.2f}")

ks = list(R)
if len(ks) == 2:
    c, p = R[ks[0]], R[ks[1]]
    print(f"\n  MODEL-MEASURED PERIOD OFFSET (presentday minus control):")
    print(f"    Siberian JJA T2m : {p['sib']-c['sib']:+.2f} K"
          f"   (HadCRUT5-chain estimate was +1.12 K)")
    print(f"    SO SW CRE        : {p['so_cre']-c['so_cre']:+.2f} W/m2"
          f"   (ERA5-TOA estimate was -1.43)")
    print(f"    global net TOA   : {p['toa']-c['toa']:+.2f} W/m2")
    print(f"\n  => The genuine boreal cold bias, period-clean, is "
          f"{p['sib']-E5:+.2f} K, not {c['sib']-E5:+.2f} K.")
