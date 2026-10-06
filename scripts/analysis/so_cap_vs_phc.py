"""Southern Ocean surface cap in the model against the PHC3 monthly climatology.

Same regions and layers as so_destabilisation.py (NEW = convects in PICAL_crunveg_kpp_gm1000 but not in
the TKE-only line; OLD = convects in both, mostly Weddell). For March and September: salinity and
temperature of 0-50 m, 50-150 m and 200-500 m, and the density step sigma0(200-300 m) - sigma0(0-50 m)
with its salinity and temperature parts. Model: PICAL_crunveg_tke_albsn082, mean of 2200-2209.
PHC3 (Steele et al. 2001, 1 degree, monthly to 1500 m) is sampled at the nearest grid point of every
region node. Caveat: winter observations under Antarctic sea ice are sparse and the climatology is
smooth, so its September values lean on few profiles.

Usage: so_cap_vs_phc.py    (reval environment; writes data/clim/so_cap_vs_phc.txt)
"""
import os
import numpy as np
import xarray as xr
import gsw

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
RUN = f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082/outdata/fesom'
ARM = f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_kpp_gm1000/outdata/fesom'
PHC = f'{P}/input/fesom2/hydrography_recom/phc3.0_monthly.nc'
YEARS = range(2200, 2210)
md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lon, lat = md['lon'].values, md['lat'].values
zi = np.abs(md['nz'].values); zm, dz = 0.5 * (zi[1:] + zi[:-1]), np.diff(zi)
area = md['nod_area'].values[0]
od = (md['ulevels_nod2D'].values == 1) & (np.abs(md['zbar_n_bottom'].values) > 2000) & (lat < -55)
mld = lambda d: np.abs(xr.open_dataset(f'{d}/MLD2.fesom.2177.nc')['MLD2'].isel(time=8).values)
ma, mt = mld(ARM), mld(RUN)
REG = {'NEW': od & (ma > 1000) & (mt < 500), 'OLD': od & (ma > 1000) & (mt > 1000)}
LAY = [('0-50', 0, 50), ('50-150', 50, 150), ('200-300', 200, 300), ('200-500', 200, 500)]


def step(s_top, t_top, s_bot, t_bot):
    sm, tm = 0.5 * (s_top + s_bot), 0.5 * (t_top + t_bot)
    rho = 1000 + gsw.sigma0(sm, tm)
    return (float(gsw.sigma0(s_bot, t_bot) - gsw.sigma0(s_top, t_top)), float(rho * gsw.beta(sm, tm, 150) * (s_bot - s_top)),
            float(-rho * gsw.alpha(sm, tm, 150) * (t_bot - t_top)))


def model(idx, w, mo):
    acc = {}
    for y in YEARS:
        T = xr.open_dataset(f'{RUN}/temp.fesom.{y}.nc')['temp'].isel(time=mo).values[idx]
        S = xr.open_dataset(f'{RUN}/salt.fesom.{y}.nc')['salt'].isel(time=mo).values[idx]
        T = np.where(S > 1, T, np.nan); S = np.where(S > 1, S, np.nan)
        for k, a, b in LAY:
            m = (zm >= a) & (zm < b)
            for nm, x in (('S', S), ('T', T)):
                v = np.nansum(x[:, m] * dz[m], 1) / np.nansum(np.where(np.isfinite(x[:, m]), dz[m], 0), 1)
                acc.setdefault(nm + k, []).append(float(np.nansum(v * w)))
    return {k: float(np.mean(v)) for k, v in acc.items()}


ph = xr.open_dataset(PHC, decode_times=False)
pz = ph['depth'].values
pdz = np.gradient(pz)


def phc(idx, w, mo):
    la = xr.DataArray(lat[idx], dims='n'); lo = xr.DataArray(lon[idx] % 360, dims='n')
    out = {}
    T = ph['temp'].isel(time=mo).sel(lat=la, lon=lo, method='nearest').values      # (depth, n)
    S = ph['salt'].isel(time=mo).sel(lat=la, lon=lo, method='nearest').values
    for k, a, b in LAY:
        m = (pz >= a) & (pz <= b)
        for nm, x in (('S', S), ('T', T)):
            v = np.nansum(x[m] * pdz[m, None], 0) / np.nansum(np.where(np.isfinite(x[m]), pdz[m, None], 0), 0)
            ok = np.isfinite(v)
            out[nm + k] = float(np.sum(v[ok] * w[ok]) / np.sum(w[ok]))
    return out


out = ['Southern Ocean surface cap: PICAL_crunveg_tke_albsn082 (2200-2209) against PHC3 monthly climatology.',
       'step = sigma0(200-300 m) - sigma0(0-50 m) [kg/m3], with its salinity and temperature parts.']
for name, reg in REG.items():
    idx = np.where(reg)[0]; w = area[idx] / area[idx].sum()
    out += ['', f'Region {name} ({area[idx].sum() / 1e12:.2f}e6 km2)',
            f"{'':10s}{'':6s}{'S 0-50':>9s}{'S 50-150':>9s}{'S 200-500':>10s}{'T 0-50':>8s}{'T 50-150':>9s}{'T 200-500':>10s}{'step':>8s}{'salt':>8s}{'temp':>8s}"]
    for mo, mn in ((2, 'March'), (8, 'September')):
        for src, f in (('model', model), ('PHC3', phc)):
            r = f(idx, w, mo)
            s = step(r['S0-50'], r['T0-50'], r['S200-300'], r['T200-300'])
            out.append(f"{mn:10s}{src:6s}{r['S0-50']:9.3f}{r['S50-150']:9.3f}{r['S200-500']:10.3f}{r['T0-50']:8.2f}{r['T50-150']:9.2f}{r['T200-500']:10.2f}{s[0]:8.3f}{s[1]:8.3f}{s[2]:8.3f}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/so_cap_vs_phc.txt', 'w').write(txt + '\n')
