"""Is the boreal-summer surface shortwave deficit a CLOUD problem, and if so what kind?

The AMIP baseline is 2.16 K too cold over Siberian land in JJA, with surface net SW
33 W/m2 below CRUNCEP3. CRUNCEP3 is itself a reanalysis-derived product, so this
re-tests the deficit against CERES EBAF Ed4.1 (satellite) and decomposes it:

    d(SWnet_all)  =  d(SWnet_clear)  +  d(SW CRE at surface)
                     ^clear-sky error   ^cloud error
                     (aerosol, water     (amount x opacity)
                      vapour, albedo)

and then splits the cloud part by comparing cloud AREA directly (model tcc vs CERES
cldarea). Area right + CRE too negative => clouds too optically thick. Area too high
=> too much cloud.

Using NET fluxes on both sides avoids any assumption about the surface albedo or the
direct/diffuse split, since CERES supplies all-sky and clear-sky net SW directly.

Land-specificity is tested by repeating over ocean in the same latitude band: if the
error is there too, it is a general cloud problem and the safe knobs are different.

CAVEAT: CERES is 2000s; the model is ~PI (1850 aerosol, 1870s SST). For CLEAR-SKY SW
the PI-vs-present aerosol difference works AGAINST a "model too dark" finding (less PI
aerosol => PI clear-sky SW should be HIGHER than CERES), so a clear-sky deficit here
would be a conservative result. Cloud amount and CRE are far less sensitive to that.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

ACC = 3600.0
REPO = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/'
AMIP = '/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_base/outdata/oifs'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
Y0, Y1 = 1872, 1879
FLUX = {'ssr', 'ssrc', 'ssrd', 'tsr', 'tsrc'}
VARS = ['ssr', 'ssrc', 'ssrd', 'tsr', 'tsrc', 'tcc', 'lcc', 'mcc', 'hcc']

# (name, lat0, lat1, lon0, lon1, surface) with lon in -180..180
REGIONS = [('Siberia land',    55, 75,   60, 180, 'land'),
           ('E. Siberia land', 55, 75,   90, 160, 'land'),
           ('NH 55-75N ocean', 55, 75, -180, 180, 'ocean')]


def load_model():
    out = {}
    for v in VARS:
        acc = []
        for y in range(Y0, Y1 + 1):
            f = f'{AMIP}/atm_remapped_1m_{v}_1m_{y}-{y}.nc'
            if not os.path.exists(f):
                continue
            ds = xr.open_dataset(f)
            a = ds[v].values                       # (12, lat, lon)
            if v in FLUX:
                a = a / ACC
            acc.append(a)
            lat, lon = ds[v].lat.values, ds[v].lon.values
            ds.close()
        out[v] = np.mean(acc, axis=0)
    return out, lat, lon


def load_ceres(lat, lon):
    """CERES monthly climatology, bilinearly interpolated to the model grid."""
    ds = xr.open_dataset(CERESF)
    # CERES lon is 0.5..359.5, so target lons below 0.5 or above 359.5 would
    # interpolate to NaN and poison any all-longitude average. Wrap the field
    # periodically first (append lon-360 at the start and lon+360 at the end).
    ds = xr.concat([ds.isel(lon=[-1]).assign_coords(lon=ds.lon.values[-1:] - 360.0),
                    ds,
                    ds.isel(lon=[0]).assign_coords(lon=ds.lon.values[:1] + 360.0)],
                   dim='lon')
    tgt_lon = np.where(lon < 0, lon + 360, lon)
    # The model grid runs to +-89.53 while CERES stops at +-89.5, so the two polar
    # rows would extrapolate to NaN. Clamp them; the 0.03 deg shift is immaterial
    # and those rows are outside every region used here.
    tgt_lat = np.clip(lat, ds.lat.values.min(), ds.lat.values.max())
    keep = ['sfc_net_sw_all_clim', 'sfc_net_sw_clr_t_clim', 'sfc_sw_down_all_clim',
            'cldarea_total_daynight_clim', 'cldtau_total_day_clim',
            'toa_cre_sw_clim', 'sfc_cre_net_sw_clim']
    out = {}
    for v in keep:
        da = ds[v].interp(lat=xr.DataArray(tgt_lat, dims='y'),
                          lon=xr.DataArray(tgt_lon, dims='x'))
        out[v] = da.values                          # (ctime, y, x)
        # cldtau is a DAYTIME retrieval, so it is genuinely absent in polar night.
        if v != 'cldtau_total_day_clim':
            assert np.isfinite(out[v]).all(), f'NaN remains in interpolated {v}'
    ds.close()
    return out


def boxmean(field, lat, lon, r, lsm, months):
    """field (12, lat, lon) -> scalar over region r and the given month indices."""
    _, la0, la1, lo0, lo1, sfc = r
    yi = (lat >= la0) & (lat <= la1)
    lon180 = ((lon + 180) % 360) - 180
    xi = (lon180 >= lo0) & (lon180 <= lo1)
    sub = field[np.ix_(months, np.where(yi)[0], np.where(xi)[0])]
    L = lsm[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    m = (L > 0.5) if sfc == 'land' else (L <= 0.5)
    w = np.cos(np.deg2rad(lat[yi]))[:, None] * np.ones(m.shape)
    w = np.where(m, w, 0.0)
    W = np.broadcast_to(w, sub.shape).copy()
    W[~np.isfinite(sub)] = 0.0          # NaN-aware (cldtau in polar night)
    if W.sum() == 0:
        return np.nan
    return float(np.nansum(sub * W) / W.sum())


mod, lat, lon = load_model()
cer = load_ceres(lat, lon)
lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values
JJA = [5, 6, 7]

print(f"AMIP {Y0}-{Y1} vs CERES EBAF Ed4.1, JJA, area-weighted\n")
hdr = f"{'':34s}" + "".join(f"{r[0]:>19s}" for r in REGIONS)
print(hdr)


def row(label, mfield, cfield, fmt="{:7.1f}"):
    cells = []
    for r in REGIONS:
        m = boxmean(mfield, lat, lon, r, lsm, JJA)
        c = boxmean(cfield, lat, lon, r, lsm, JJA) if cfield is not None else np.nan
        cells.append(f"{fmt.format(m)}/{fmt.format(c)}/{m-c:+6.1f}" if cfield is not None
                     else f"{fmt.format(m):>19s}")
    print(f"{label:34s}" + "".join(f"{c:>19s}" for c in cells))


print("  (model / CERES / difference)\n")
row('SW net at surface, all-sky', mod['ssr'], cer['sfc_net_sw_all_clim'])
row('SW net at surface, clear-sky', mod['ssrc'], cer['sfc_net_sw_clr_t_clim'])
row('SW CRE at surface', mod['ssr'] - mod['ssrc'], cer['sfc_cre_net_sw_clim'])
row('SW down at surface, all-sky', mod['ssrd'], cer['sfc_sw_down_all_clim'])
row('TOA SW CRE', mod['tsr'] - mod['tsrc'], cer['toa_cre_sw_clim'])
row('cloud area [%]', mod['tcc'] * 100, cer['cldarea_total_daynight_clim'])
print()
for r in REGIONS[:1]:
    pass
print(f"{'model low/mid/high cloud [%]':34s}" + "".join(
    f"{boxmean(mod['lcc'],lat,lon,r,lsm,JJA)*100:5.1f}/"
    f"{boxmean(mod['mcc'],lat,lon,r,lsm,JJA)*100:5.1f}/"
    f"{boxmean(mod['hcc'],lat,lon,r,lsm,JJA)*100:5.1f}".rjust(19) for r in REGIONS))
print(f"{'CERES cloud optical depth':34s}" + "".join(
    f"{boxmean(cer['cldtau_total_day_clim'],lat,lon,r,lsm,JJA):19.1f}" for r in REGIONS))

print("\n\nMonthly, Siberia land (model / CERES / diff):")
print(f"  {'month':>5s} {'SWnet all':>22s} {'SWnet clear':>22s} {'SW CRE sfc':>22s} {'cloud area %':>22s}")
for mi in range(12):
    r = REGIONS[0]
    f = lambda M, C: (boxmean(M, lat, lon, r, lsm, [mi]), boxmean(C, lat, lon, r, lsm, [mi]))
    a = f(mod['ssr'], cer['sfc_net_sw_all_clim'])
    b = f(mod['ssrc'], cer['sfc_net_sw_clr_t_clim'])
    c = f(mod['ssr'] - mod['ssrc'], cer['sfc_cre_net_sw_clim'])
    d = f(mod['tcc'] * 100, cer['cldarea_total_daynight_clim'])
    print(f"  {mi+1:5d} " + " ".join(f"{x[0]:7.1f}/{x[1]:6.1f}/{x[0]-x[1]:+6.1f}" for x in (a, b, c, d)))
