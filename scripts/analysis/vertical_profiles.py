"""Where in the column do the biases live? Model minus ERA5 vertical profiles.

Every metric the campaign has used is at the surface or TOA, which says a bias exists but
not where it originates. This asks:

  * Is the boreal JJA cold bias confined to the boundary layer -- a surface-exchange problem,
    which is what the F-series (z0h, LAI, cover, stomatal resistance) assumes -- or does it
    extend through the troposphere, which would make it a circulation or radiation problem
    and the F-series largely beside the point?
  * Is the Southern Ocean moist/cloud bias in LOW cloud? D2b's pressure gate (reduce INP only
    below 700 hPa) beat the ungated D2a on SO SW RMSE, which implies the SO signal is a
    low-cloud phenomenon. RH profiles test that independently.
  * Do the tropics show a deep moist/warm bias consistent with their 2.5 W/m2 TOA deficit?

Model does not write cloud on pressure levels, so RH is the cloud proxy throughout.

WARNING -- BROKEN MODEL DIAGNOSTIC. The model's pressure-level relative humidity field `r`
is IDENTICALLY ZERO in every run checked (0 nonzero of 17,510,400 values in both
amip_presentday 1995 and amip_pi_base 1900), while `q` on the same levels is fine. The XIOS
`r` output is therefore dead and must not be used. RH here is computed from t and q with
Bolton's saturation formula over water, applied IDENTICALLY to model and ERA5, so any
convention difference cancels in the difference.

PRIMARY comparison is amip_presentday (1990-2014) vs ERA5 (1990-2014): same epoch, so no
reference-period offset. The PI comparison (control 1872-1915 vs ERA5 1940-1969) is shown
alongside but carries the residual epoch offset that amip_presentday measured at +0.42 K on
Siberian JJA T2m -- so read its absolute values with that caveat, and prefer the SHAPE of the
profile to its offset.

Regions use the same boxes and the same land/sea mask as eval_round10_A.py, applied to both
datasets after ERA5 has been regridded to the model grid, so the masking is identical.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

W = '/tmp/vprof'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values

REGIONS = {
    'Siberia land  JJA': dict(la=(55, 75), lo=(60, 180), sfc='land', season='JJA'),
    'SO 45-65S ocn DJF': dict(la=(-65, -45), lo=(-180, 180), sfc='ocean', season='DJF'),
    'SO 45-65S ocn ANN': dict(la=(-65, -45), lo=(-180, 180), sfc='ocean', season='ANN'),
    'Tropics 30S-30N  ': dict(la=(-30, 30), lo=(-180, 180), sfc='all', season='ANN'),
    'Global           ': dict(la=(-90, 90), lo=(-180, 180), sfc='all', season='ANN'),
}
VARS = [('t', 'T [K]'), ('rh', 'RH [%]'), ('q', 'q [g/kg]')]
SHOW = [100000, 92500, 85000, 70000, 50000, 30000, 20000, 10000]   # Pa, surface -> 100 hPa


def load(kind, tag, v, season):
    f = f'{W}/{kind}_{tag}_{v}_{season}.nc'
    if not os.path.exists(f):
        return None, None
    d = xr.open_dataset(f)
    name = [k for k in d.data_vars
            if 'bnds' not in k and 'bounds' not in k][0]
    a = d[name]
    lev = [c for c in a.dims if 'lev' in c or 'plev' in c][0]
    arr = a.squeeze().transpose(lev, ...).values
    levs = d[lev].values.astype(float)
    lat = d['lat'].values
    lon = d['lon'].values
    d.close()
    return (arr, levs, lat, lon), name


def rh_from_tq(t, q, levs):
    """RH [%] from T [K], q [kg/kg] and pressure [Pa]; Bolton (1980) over water.

    Applied identically to model and ERA5, so the choice of saturation convention
    (water-only here, versus the IFS mixed-phase definition) cancels in the difference.
    """
    p = levs[:, None, None]
    es = 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))
    qs = 0.622 * es / np.maximum(p - 0.378 * es, 1.0)
    return 100.0 * q / np.maximum(qs, 1e-12)


def load_field(kind, tag, v, season):
    """Return (field, levs, lat, lon) for a real or derived variable."""
    if v != 'rh':
        r = load(kind, tag, v, season)[0]
        return r
    T = load(kind, tag, 't', season)[0]
    Q = load(kind, tag, 'q', season)[0]
    if T is None or Q is None:
        return None
    t, levs, lat, lon = T
    q = Q[0]
    return rh_from_tq(t, q, levs), levs, lat, lon


def profile(arr, levs, lat, lon, la, lo, sfc):
    yi = (lat >= la[0]) & (lat <= la[1])
    l180 = ((lon + 180) % 360) - 180
    xi = (l180 >= lo[0]) & (l180 <= lo[1])
    L = lsm[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    m = np.ones(L.shape, bool) if sfc == 'all' else (L > 0.5 if sfc == 'land' else L <= 0.5)
    w = np.broadcast_to(np.cos(np.deg2rad(lat[yi]))[:, None], L.shape)
    out = {}
    for k, p in enumerate(levs):
        sub = arr[k][np.ix_(np.where(yi)[0], np.where(xi)[0])]
        ok = m & np.isfinite(sub)
        out[int(round(p))] = np.average(sub[ok], weights=w[ok]) if ok.any() else np.nan
    return out


for tag, title in (('pd', 'PRESENT DAY  model 1990-2014 vs ERA5 1990-2014  (period-clean)'),
                   ('pi', 'PI EPOCH     model 1872-1915 vs ERA5 1940-1969  (epoch offset applies)')):
    print(f'\n{"="*100}\n{title}\n{"="*100}')
    for rname, R in REGIONS.items():
        rows = []
        for v, vlabel in VARS:
            M = load_field('model', tag, v, R['season'])
            E = load_field('era5', tag, v, R['season'])
            if M is None or E is None:
                continue
            pm = profile(*M, R['la'], R['lo'], R['sfc'])
            pe = profile(*E, R['la'], R['lo'], R['sfc'])
            scale = 1000.0 if v == 'q' else 1.0     # kg/kg -> g/kg
            rows.append((vlabel, {p: (pm[p] - pe[p]) * scale for p in pm if p in pe}))
        if not rows:
            continue
        print(f'\n  {rname}   model - ERA5')
        print('    ' + 'level [hPa]'.ljust(14) + ''.join(f'{p//100:>9d}' for p in SHOW))
        for vlabel, d in rows:
            print('    ' + vlabel.ljust(14) + ''.join(
                f'{d[p]:9.2f}' if p in d and np.isfinite(d[p]) else f'{"--":>9s}' for p in SHOW))
