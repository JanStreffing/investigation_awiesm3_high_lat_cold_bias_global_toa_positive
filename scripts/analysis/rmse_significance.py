"""Are the deep-water SW RMSE differences real, or the same noise that ate the boreal table?

noise_floor.py tested regional MEANS and found the energy target well resolved and the
boreal target not resolved at all. It did not test spatial RMSE, which is a different
statistic -- and the deep-water SW RMSEs are the metric the campaign calls "the priority",
because SW errors where deep water forms set the coupled ocean's initial state.

One claim in particular is load-bearing and unchecked: that combining A1b with B2 CANCELS
B2's subpolar-Atlantic damage (+1.226 alone -> +0.071 combined). If that is noise, the case
that the energy and boreal targets do not trade against each other in that basin collapses.

Method: same run x year ANOVA as noise_floor.py. For each run and each year, compute the
area-weighted spatial RMSE of that YEAR's mean surface-SW field against the CERES
climatology, then decompose

    X[run, year] = mu + a[run] + g[year] + eps[run, year]

and report SE of a run-minus-control difference = sd(eps) * sqrt(2/n_years).

NOTE on interpretation: the headline table quotes RMSE of the 4-year MEAN field, which is
smaller than the mean of per-year RMSEs because averaging suppresses internal variability.
The per-year RMSE used here is the right replicate for significance testing -- it asks
whether the run effect on RMSE exceeds year-to-year noise -- but its absolute values sit
above the headline numbers by construction. Compare differences, not absolute levels.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

ACC = 3600.0
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
OBS = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
Y0, Y1 = 1872, 1915
YEARS = list(range(Y0, Y1 + 1))

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
        ('ABB8 A1b+B2+B8', 'amip_ABB8'),
        ('C1 rlam75', 'amip_C1_rlam75'),
        ('C2 rlam40', 'amip_C2_rlam40'),
        ('E1 lamsk2.5', 'amip_E1_lamsk25')]

BOXES = {'SO 45-65S': (-65, -45, -180, 180), 'subpolar N Atl': (50, 65, -60, 0),
         'Nordic Seas': (65, 80, -20, 20), 'global': None}

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def ceres(lat, lon):
    ds = xr.open_dataset(OBS)
    ds = xr.concat([ds.isel(lon=[-1]).assign_coords(lon=ds.lon.values[-1:] - 360.0), ds,
                    ds.isel(lon=[0]).assign_coords(lon=ds.lon.values[:1] + 360.0)], dim='lon')
    tl = np.clip(lat, -89.5, 89.5)
    tlo = np.where(lon < 0, lon + 360, lon)
    return ds['sfc_net_sw_all_clim'].interp(lat=xr.DataArray(tl, dims='y'),
                                            lon=xr.DataArray(tlo, dims='x')).values


def rmse(model2d, obs2d, lat, lon, box, ocean_only):
    w = np.cos(np.deg2rad(lat))[:, None] * np.ones(model2d.shape)
    if box is not None:
        la0, la1, lo0, lo1 = box
        l180 = ((lon + 180) % 360) - 180
        m = np.zeros(model2d.shape, bool)
        m[np.ix_((lat >= la0) & (lat <= la1), (l180 >= lo0) & (l180 <= lo1))] = True
        if ocean_only:
            m &= (lsm <= 0.5)
        w = np.where(m, w, 0.0)
    ok = np.isfinite(model2d) & np.isfinite(obs2d)
    w = np.where(ok, w, 0.0)
    return float(np.sqrt((w * (model2d - obs2d) ** 2).sum() / w.sum()))


def per_year_rmse(run, C, lat, lon):
    out = {k: [] for k in BOXES}
    for y in YEARS:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_ssr_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            return None
        ds = xr.open_dataset(f)
        fld = ds['ssr'].values.mean(axis=0) / ACC
        ds.close()
        for k, b in BOXES.items():
            out[k].append(rmse(fld, C, lat, lon, b, b is not None))
    return out


# grid + CERES from the first available run
lat = lon = C = None
for _, run in RUNS:
    f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_ssr_1m_{Y0}-{Y0}.nc'
    if os.path.exists(f):
        ds = xr.open_dataset(f); lat, lon = ds['ssr'].lat.values, ds['ssr'].lon.values; ds.close()
        C = ceres(lat, lon); break

labs, data = [], {k: [] for k in BOXES}
for lab, run in RUNS:
    r = per_year_rmse(run, C, lat, lon)
    if r is None:
        print(f'  !! {lab}: incomplete, skipped'); continue
    labs.append(lab)
    for k in BOXES:
        data[k].append(r[k])

ctl = labs.index('control')
print(f'\nPer-year surface-SW RMSE vs CERES, {len(labs)} runs x {len(YEARS)} yr '
      f'({Y0}-{Y1}). Run x year ANOVA.\n')
for k in BOXES:
    X = np.array(data[k]); nr, ny = X.shape
    mu = X.mean(); a = X.mean(1) - mu; g = X.mean(0) - mu
    eps = X - (mu + a[:, None] + g[None, :])
    sd = np.sqrt((eps ** 2).sum() / ((nr - 1) * (ny - 1)))
    sed = sd * np.sqrt(2.0 / ny)
    print(f'{k:16s}  sd(eps)={sd:6.3f}  SE_diff={sed:6.3f}  95%thr=±{1.96*sed:5.3f}')
    sig = [(labs[i], X[i].mean() - X[ctl].mean()) for i in range(nr)
           if i != ctl and abs(X[i].mean() - X[ctl].mean()) > 1.96 * sed]
    print('   significant: ' + (', '.join(f'{l}{v:+.2f}' for l, v in sig) if sig else 'NONE'))
    # the load-bearing claim
    if k == 'subpolar N Atl':
        for nm in ('B2 convi=25', 'AB ovl+convi', 'ABB8 A1b+B2+B8'):
            if nm in labs:
                d = X[labs.index(nm)].mean() - X[ctl].mean()
                print(f'   -> {nm:16s} {d:+.3f}  t={d/sed:+.2f}  '
                      f'{"SIGNIFICANT" if abs(d/sed) > 1.96 else "within noise"}')
    print()
