"""Round 10 prong A evaluation: A1 (Southern Ocean) and A2 (boreal land) vs control.

Runs (all AMIP TCO95, 1850 GHG, observed SST 1870s; evaluated on 1872-1875):
  control  amip_pi_base            RCL_OVERLAPLIQICE=0.65, RCL_KK_CLOUD_NUM_LAND=300
  A1a      amip_A1_overlap01       RCL_OVERLAPLIQICE=0.10   (EC-Earth4's tuned value)
  A1b      amip_A1_overlap035      RCL_OVERLAPLIQICE=0.35
  A2       amip_A2_kknumland150    RCL_KK_CLOUD_NUM_LAND=150 (source change)
  expA     amip_expA_rvrsmin500    RVRSMIN(3,4)=500          (partitioning test)

Two falsifiable predictions being tested:
  * A1 should leave BOREAL JJA essentially unchanged -- the WBF deposition term only
    acts where supercooled liquid coexists with ice (RTHOMO <= T < RTT-5 K), and
    boreal-land summer BL cloud is warm (>268 K) and never enters that window.
  * A2 should leave the SOUTHERN OCEAN unchanged -- RCL_KK_CLOUD_NUM_LAND is used
    only in the PLSM>0.5 branch.
If either bleeds, the design's separability argument is wrong.

Guardrails reported for every run: global net TOA and surface flux (with snow
enthalpy), tropics, and NH-SH albedo asymmetry (which the literature warns A1 and A2
push the SAME way).
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

ACC = 3600.0
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
OBS = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
Y0, Y1 = 1872, 1875
FLUX = {'tsr', 'ttr', 'tsrc', 'ttrc', 'ssr', 'str', 'sshf', 'slhf', 'sf', 'ssrd', 'tisr'}
VARS = ['tsr', 'ttr', 'tsrc', 'ttrc', 'tisr', 'ssr', 'str', 'sshf', 'slhf', 'sf', 'tcc', '2t']

RUNS = [('control', 'amip_pi_base'), ('A1a ovl=0.10', 'amip_A1_overlap01'),
        ('A1b ovl=0.35', 'amip_A1_overlap035'), ('A2 KKland=150', 'amip_A2_kknumland150'),
        ('expA rvrs=500', 'amip_expA_rvrsmin500'),
        ('A1c depth1500', 'amip_A1c_depliqdepth1500'),
        ('B1 detrpen.45', 'amip_B1_detrpen045'),
        ('B2 convi=25', 'amip_B2_clddiffconvi25'),
        ('AB ovl+convi', 'amip_AB_ovl035_convi25'),
        ('B3 clddiff', 'amip_B3_clddiff15e6'),
        ('B4 entshalp3', 'amip_B4_entshalp3'),
        ('B5 capdcycl0', 'amip_B5_capdcycl0'),
        ('B6 lcritsnow', 'amip_B6_lcritsnow1e5'),
        ('B7 rvice.22', 'amip_B7_rvice022'),
        ('B8 lamsk5', 'amip_B8_lamsk5'),
        # --- started 2026-07-29. C1/C2/E1 sit on the B8 binary, so their
        # reference is 'B8 lamsk5', not 'control'.
        ('ABB8 A1b+B2+B8', 'amip_ABB8'),
        ('C1 rlam75', 'amip_C1_rlam75'),
        ('C2 rlam40', 'amip_C2_rlam40'),
        ('E1 lamsk2.5', 'amip_E1_lamsk25')]

ERA5 = '/work/ab0246/a270092/obs/era5/netcdf/T2M.nc'


def load(run):
    out = {}
    for v in VARS:
        acc = []
        for y in range(Y0, Y1 + 1):
            f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc'
            if not os.path.exists(f):
                return None, None, None
            ds = xr.open_dataset(f)
            a = ds[v].values
            if v in FLUX:
                a = a / ACC
            acc.append(a)
            lat, lon = ds[v].lat.values, ds[v].lon.values
            ds.close()
        out[v] = np.mean(acc, axis=0)
    return out, lat, lon


def ceres(lat, lon):
    ds = xr.open_dataset(OBS)
    ds = xr.concat([ds.isel(lon=[-1]).assign_coords(lon=ds.lon.values[-1:] - 360.0), ds,
                    ds.isel(lon=[0]).assign_coords(lon=ds.lon.values[:1] + 360.0)], dim='lon')
    tl = np.clip(lat, -89.5, 89.5)
    tlo = np.where(lon < 0, lon + 360, lon)
    keep = ['toa_cre_sw_clim', 'cldarea_total_daynight_clim', 'toa_net_all_clim',
            'sfc_net_sw_all_clim', 'sfc_cre_net_sw_clim']
    return {v: ds[v].interp(lat=xr.DataArray(tl, dims='y'),
                            lon=xr.DataArray(tlo, dims='x')).values for v in keep}


def sel(f, lat, lon, lsm, la, lo=(-180, 180), sfc='all', months=None):
    yi = (lat >= la[0]) & (lat <= la[1])
    l180 = ((lon + 180) % 360) - 180
    xi = (l180 >= lo[0]) & (l180 <= lo[1])
    mo = months if months is not None else list(range(12))
    sub = f[np.ix_(mo, np.where(yi)[0], np.where(xi)[0])]
    L = lsm[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    m = np.ones(L.shape, bool) if sfc == 'all' else (L > 0.5 if sfc == 'land' else L <= 0.5)
    w = np.cos(np.deg2rad(lat[yi]))[:, None] * np.ones(m.shape)
    w = np.where(m, w, 0.0)
    return float((sub * w).sum() / (w.sum() * len(mo)))


def era5_t2m(lat, lon):
    """ERA5 T2m annual climatology on the model grid. NOTE the model is ~PI
    (1850 GHG, 1870s SST) and ERA5 is present-day, so an expected cold offset of
    order 0.5-1 K is folded into this RMSE -- use it to RANK runs, not as an
    absolute skill score."""
    ds = xr.open_dataset(ERA5)
    v = 't2m' if 't2m' in ds else list(ds.data_vars)[0]
    da = ds[v].mean('time') if 'time' in ds[v].dims else ds[v]
    la = 'latitude' if 'latitude' in da.dims else 'lat'
    lo = 'longitude' if 'longitude' in da.dims else 'lon'
    da = da.rename({la: 'lat', lo: 'lon'}).sortby('lat')
    da = xr.concat([da.isel(lon=[-1]).assign_coords(lon=da.lon.values[-1:] - 360.0), da,
                    da.isel(lon=[0]).assign_coords(lon=da.lon.values[:1] + 360.0)], dim='lon')
    out = da.interp(lat=xr.DataArray(np.clip(lat, da.lat.values.min(), da.lat.values.max()), dims='y'),
                    lon=xr.DataArray(np.where(lon < 0, lon + 360, lon), dims='x')).values
    ds.close()
    return out


def rmse(model2d, obs2d, lat, lon=None, box=None, lsm=None, sfc='all'):
    """Area-weighted spatial RMSE of an annual-mean field, optionally over a box.

    The regional versions matter more than the global one: SW biases where DEEP WATER
    FORMS set the coupled ocean's initial state, and errors there produce long,
    unpredictable coupled spin-up drift. Boxes used:
      SO     45-65S ocean            -- AABW/AAIW formation
      SPNA   50-65N, 60W-0, ocean    -- NADW, subpolar North Atlantic
      NORDIC 65-80N, 20W-20E, ocean  -- Nordic Seas overflow
    """
    w = np.cos(np.deg2rad(lat))[:, None] * np.ones(model2d.shape)
    if box is not None:
        la0, la1, lo0, lo1 = box
        l180 = ((lon + 180) % 360) - 180
        m = np.zeros(model2d.shape, bool)
        m[np.ix_((lat >= la0) & (lat <= la1), (l180 >= lo0) & (l180 <= lo1))] = True
        if lsm is not None and sfc != 'all':
            m &= (lsm > 0.5) if sfc == 'land' else (lsm <= 0.5)
        w = np.where(m, w, 0.0)
    ok = np.isfinite(model2d) & np.isfinite(obs2d)
    w = np.where(ok, w, 0.0)
    return float(np.sqrt((w * (model2d - obs2d) ** 2).sum() / w.sum()))


lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values
JJA = [5, 6, 7]
E5 = None
res, C = {}, None
for lab, run in RUNS:
    d, _lat, _lon = load(run)
    if d is None:
        print(f'  !! {lab}: output incomplete, skipped'); continue
    lat, lon = _lat, _lon      # keep the last GOOD grid; a failed load must not clobber it
    if C is None:
        C = ceres(lat, lon)
        E5 = era5_t2m(lat, lon)
    net = d['tsr'] + d['ttr']
    swcre = d['tsr'] - d['tsrc']
    snow = d['sf'] * 333_550_000.0
    sfc = d['ssr'] + d['str'] + d['sshf'] + d['slhf'] - snow
    alb_n = 1 - sel(d['tsr'], lat, lon, lsm, (0, 90)) / sel(d['tisr'], lat, lon, lsm, (0, 90))
    alb_s = 1 - sel(d['tsr'], lat, lon, lsm, (-90, 0)) / sel(d['tisr'], lat, lon, lsm, (-90, 0))
    res[lab] = dict(
        so_cre=sel(swcre, lat, lon, lsm, (-65, -45), sfc='ocean'),
        so_cld=sel(d['tcc'] * 100, lat, lon, lsm, (-65, -45), sfc='ocean'),
        so_net=sel(net, lat, lon, lsm, (-65, -45), sfc='ocean'),
        sib_t=sel(d['2t'], lat, lon, lsm, (55, 75), (60, 180), 'land', JJA) - 273.15,
        sib_sw=sel(d['ssr'], lat, lon, lsm, (55, 75), (60, 180), 'land', JJA),
        sib_cld=sel(d['tcc'] * 100, lat, lon, lsm, (55, 75), (60, 180), 'land', JJA),
        g_toa=sel(net, lat, lon, lsm, (-90, 90)),
        g_sfc=sel(sfc, lat, lon, lsm, (-90, 90)),
        trop=sel(net, lat, lon, lsm, (-30, 30)),
        dalb=alb_n - alb_s,
        rmse_sw=rmse(d['ssr'].mean(0), C['sfc_net_sw_all_clim'].mean(0), lat),
        rmse_so=rmse(d['ssr'].mean(0), C['sfc_net_sw_all_clim'].mean(0), lat, lon,
                     (-65, -45, -180, 180), lsm, 'ocean'),
        rmse_spna=rmse(d['ssr'].mean(0), C['sfc_net_sw_all_clim'].mean(0), lat, lon,
                       (50, 65, -60, 0), lsm, 'ocean'),
        rmse_nordic=rmse(d['ssr'].mean(0), C['sfc_net_sw_all_clim'].mean(0), lat, lon,
                         (65, 80, -20, 20), lsm, 'ocean'),
        rmse_t2m=rmse(d['2t'].mean(0), E5, lat))
obs = dict(so_cre=sel(C['toa_cre_sw_clim'], lat, lon, lsm, (-65, -45), sfc='ocean'),
           so_cld=sel(C['cldarea_total_daynight_clim'], lat, lon, lsm, (-65, -45), sfc='ocean'),
           so_net=sel(C['toa_net_all_clim'], lat, lon, lsm, (-65, -45), sfc='ocean'),
           sib_sw=sel(C['sfc_net_sw_all_clim'], lat, lon, lsm, (55, 75), (60, 180), 'land', JJA),
           sib_cld=sel(C['cldarea_total_daynight_clim'], lat, lon, lsm, (55, 75), (60, 180), 'land', JJA),
           trop=sel(C['toa_net_all_clim'], lat, lon, lsm, (-30, 30)))

c = res['control']
rows = [('SOUTHERN OCEAN 45-65S ocean (A1 target)', None),
        ('  TOA SW CRE [W/m2]', 'so_cre'), ('  cloud area [%]', 'so_cld'), ('  net TOA [W/m2]', 'so_net'),
        ('SIBERIA land JJA (A2 target)', None),
        ('  2m temperature [C]', 'sib_t'), ('  surface net SW [W/m2]', 'sib_sw'),
        ('  cloud area [%]', 'sib_cld'),
        ('GLOBAL / guardrails', None),
        ('  net TOA [W/m2]', 'g_toa'), ('  surface flux [W/m2]', 'g_sfc'),
        ('  tropics net TOA [W/m2]', 'trop'), ('  NH-SH albedo', 'dalb'),
        ('SW RMSE vs CERES -- DEEP WATER FORMATION', None),
        ('  Southern Ocean 45-65S', 'rmse_so'), ('  subpolar N Atl 50-65N', 'rmse_spna'),
        ('  Nordic Seas 65-80N', 'rmse_nordic'),
        ('GLOBAL RMSE (rank, not skill)', None),
        ('  surface net SW vs CERES', 'rmse_sw'), ('  T2m vs ERA5 [K]', 'rmse_t2m')]
labs = [l for l, _ in RUNS if l in res]
print(f"\nAMIP {Y0}-{Y1}. Each column: value (change vs control).  CERES where available.\n")
print(f"{'':34s}" + "".join(f"{l:>19s}" for l in labs) + f"{'CERES':>10s}")
for name, k in rows:
    if k is None:
        print(f"\n{name}"); continue
    cells = []
    for l in labs:
        v = res[l][k]
        cells.append(f"{v:8.3f}" if l == 'control' else f"{v:8.3f} ({v-c[k]:+6.3f})")
    o = f"{obs[k]:10.2f}" if k in obs else f"{'--':>10s}"
    print(f"{name:34s}" + "".join(f"{x:>19s}" for x in cells) + o)
