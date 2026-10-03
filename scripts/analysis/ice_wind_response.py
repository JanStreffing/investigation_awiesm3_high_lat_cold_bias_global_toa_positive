"""Is the Antarctic winter drift wrong because of the wind, or because of how the ice answers it?

1. Pattern: JJA sea-level pressure, Amundsen Sea Low (ASL) depth, position and relative central
   pressure (ASL minimum minus the 60-70S... sector mean, as in Hosking et al. 2013), model against
   NCEP-2 (independent of IFS). Plus JJA-mean 10 m wind over each ice sector.
2. Response: complex regression of daily ice drift on daily 10 m wind per sector,
   drift = a * wind, a = |a| exp(i theta): |a| = wind factor, theta = turning angle (positive =
   ice to the left of the wind). Model: FESOM daily uice/vice against OpenIFS daily 10u/10v. Obs:
   OSI-455 satellite vectors (flags 20-22, 30) against NCEP-2 daily 10 m wind, 2005-2014.
   Both on the OSI 75 km grid, over cells with OSI vectors on >= 50 % of days; model also a_ice >= 0.15.
3. Forcing divergence: divergence of the 10 m wind on the same grid and stencil.
Usage: python ice_wind_response.py <run outdata> <y0> <y1> <label>
"""
import sys, glob, numpy as np, xarray as xr, pyproj, warnings
warnings.filterwarnings('ignore')
R, y0, y1, label = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
NC = '/albedo/work/user/jstreffi/obs/NCEP2/drift_compare'
OBS = sorted(glob.glob('/albedo/work/user/jstreffi/obs/osisaf_drift/sh/*.nc'))
GOOD = [20, 21, 22, 30]; MON = [6, 7, 8]; DX = 75e3

g0 = xr.open_dataset(OBS[0]); lat, lon = g0['lat'].values, g0['lon'].values
proj = pyproj.Proj(g0['Lambert_Azimuthal_Equal_Area'].attrs['proj4_string'])
x0, y0g = proj(lon, lat); xe, ye = proj(lon + 0.01, lat); xs, ys = proj(lon, lat - 0.01)
ex, ey = xe - x0, ye - y0g; nx, ny = x0 - xs, y0g - ys
ne, nn = np.hypot(ex, ey), np.hypot(nx, ny); ex, ey, nx, ny = ex / ne, ey / ne, nx / nn, ny / nn
sx = 1.0 if float(g0.xc[1] - g0.xc[0]) > 0 else -1.0
sy = 1.0 if float(g0.yc[1] - g0.yc[0]) < 0 else -1.0
L = np.mod(lon + 180, 360) - 180
SECT = {'all': np.ones_like(lon, bool), 'Weddell': (L >= -60) & (L < 20), 'Indian': (L >= 20) & (L < 90),
        'Pacific': (L >= 90) & (L < 160), 'Ross': (L >= 160) | (L < -130), 'Amund-Bell': (L >= -130) & (L < -60)}

def to_grid(ue, vn): return ue * ex + vn * nx, ue * ey + vn * ny
def to_en(ug, vg): return ug * ex + vg * ey, ug * nx + vg * ny
def divergence(u, v, ok):
    div = np.full(u.shape, np.nan)
    c = ok[1:-1, 1:-1] & ok[1:-1, 2:] & ok[1:-1, :-2] & ok[2:, 1:-1] & ok[:-2, 1:-1]
    div[1:-1, 1:-1] = np.where(c, sx * (u[1:-1, 2:] - u[1:-1, :-2]) / (2 * DX) + sy * (v[:-2, 1:-1] - v[2:, 1:-1]) / (2 * DX), np.nan)
    return div
def interp(field, glat, glon):
    la = np.asarray(glat, float); lo = np.mod(np.asarray(glon, float), 360)
    o = np.argsort(lo); lo = lo[o]; field = field[:, o]
    lo = np.concatenate([lo[-1:] - 360, lo, lo[:1] + 360]); field = np.concatenate([field[:, -1:], field, field[:, :1]], 1)
    if la[0] > la[-1]: la = la[::-1]; field = field[::-1]
    pl, pn = lat.ravel(), np.mod(lon.ravel(), 360)
    i = np.clip(np.searchsorted(la, pl) - 1, 0, len(la) - 2); j = np.clip(np.searchsorted(lo, pn) - 1, 0, len(lo) - 2)
    wy = (pl - la[i]) / (la[i + 1] - la[i]); wx = (pn - lo[j]) / (lo[j + 1] - lo[j])
    return (field[i, j] * (1 - wx) * (1 - wy) + field[i, j + 1] * wx * (1 - wy) + field[i + 1, j] * (1 - wx) * wy
            + field[i + 1, j + 1] * wx * wy).reshape(lat.shape)

# ---- observations: OSI drift + NCEP-2 daily wind ----
ok_all = []; OD = []; OW = []
ncu = {}; ncv = {}
for f in OBS:
    with xr.open_dataset(f) as o:
        t = o.time.values[0]; yr = int(str(t)[:4])
        if int(str(t)[5:7]) not in MON: continue
        fl = o['status_flag'].values[0]; ok = np.isin(fl, GOOD)
        dt = (o['t1'].values[0] - o['t0'].values[0]).astype('timedelta64[s]').astype(float)
        dt = np.where(np.isfinite(dt) & (dt > 0), dt, 86400.0)
        ug, vg = o['dX'].values[0] * 1e3 / dt, o['dY'].values[0] * 1e3 / dt
        ok &= np.isfinite(ug) & np.isfinite(vg)
    if yr not in ncu:
        ncu[yr] = xr.open_dataset(f'{NC}/uwnd.10m.gauss.{yr}.nc')['uwnd']; ncv[yr] = xr.open_dataset(f'{NC}/vwnd.10m.gauss.{yr}.nc')['vwnd']
    tsel = np.datetime64(str(t)[:10])     # drift ends at 12 UTC; daily NCEP mean of the start day
    uw = ncu[yr].sel(time=tsel - np.timedelta64(1, 'D'), method='nearest').squeeze().values
    vw = ncv[yr].sel(time=tsel - np.timedelta64(1, 'D'), method='nearest').squeeze().values
    glat, glon = ncu[yr]['lat'].values, ncu[yr]['lon'].values
    ue_w, vn_w = interp(uw, glat, glon), interp(vw, glat, glon)
    ue_d, vn_d = to_en(ug, vg)
    OD.append(np.where(ok, ue_d + 1j * vn_d, np.nan)); OW.append(ue_w + 1j * vn_w); ok_all.append(ok)
region = np.mean(ok_all, 0) >= 0.5
OD, OW = np.stack(OD), np.stack(OW)

# ---- model ----
fe, oi = f'{R}/fesom', f'{R}/oifs'
MD = []; MW = []; MOK = []
for y in range(y0, y1 + 1):
    U = xr.open_dataset(f'{fe}/uice.fesom.gr.{y}.nc', decode_times=False); V = xr.open_dataset(f'{fe}/vice.fesom.gr.{y}.nc', decode_times=False)
    A = xr.open_dataset(f'{fe}/a_ice.fesom.gr.{y}.nc', decode_times=False)
    W10u = xr.open_dataset(f'{oi}/atm_remapped_1d_10u_{y}-{y}.nc', use_cftime=True); W10v = xr.open_dataset(f'{oi}/atm_remapped_1d_10v_{y}-{y}.nc', use_cftime=True)
    wu, wv = W10u['10u'], W10v['10v']; wt = 'time_counter' if 'time_counter' in wu.dims else 'time'
    mon = wu[wt].dt.month.values
    for t in np.where(np.isin(mon, MON))[0]:
        if t >= U.sizes['time']: continue
        di = interp(np.nan_to_num(U['uice'].values[t]), U.lat.values, U.lon.values) + 1j * interp(np.nan_to_num(V['vice'].values[t]), U.lat.values, U.lon.values)
        ai = interp(np.nan_to_num(A['a_ice'].values[t]), U.lat.values, U.lon.values)
        w = interp(wu.values[t], wu.lat.values, wu.lon.values) + 1j * interp(wv.values[t], wu.lat.values, wu.lon.values)
        ok = region & (ai >= 0.15)
        MD.append(np.where(ok, di, np.nan)); MW.append(w); MOK.append(ok)
MD, MW, MOK = np.stack(MD), np.stack(MW), np.stack(MOK)

def fit(D, W, sel):
    m = sel[None] & np.isfinite(D) & np.isfinite(W)
    d, w = D[m], W[m]
    a = np.sum(d * np.conj(w)) / np.sum(np.abs(w) ** 2)
    r2 = 1 - np.sum(np.abs(d - a * w) ** 2) / np.sum(np.abs(d - d.mean()) ** 2)
    return abs(a) * 100, np.degrees(np.angle(a)), np.mean(np.abs(w)), np.mean(np.abs(d)) * 100, r2
print(f'{label} {y0}-{y1} JJA.  Model: FESOM drift vs OpenIFS 10 m wind.  Obs: OSI-455 2005-2014 vs NCEP-2 10 m wind')
print(f'{"sector":11s} | {"wind m/s":>8s} {"drift cm/s":>10s} {"factor %":>8s} {"turn deg":>8s} {"R2":>5s} | obs {"wind":>6s} {"drift":>6s} {"factor":>7s} {"turn":>6s} {"R2":>5s} | mean wind E/N model  obs     | wind div %/day mod obs')
for s, sel in SECT.items():
    sel = sel & region
    fm = fit(MD, MW, sel); fo = fit(OD, OW, sel)
    mwe = np.nanmean(np.where(MOK & sel[None], MW, np.nan)); owe = np.nanmean(np.where(np.isfinite(OD) & sel[None], OW, np.nan))
    wdm = np.nanmean([np.nanmean(np.where(sel, divergence(*to_grid(MW[k].real, MW[k].imag), region), np.nan)) for k in range(0, len(MW), 3)]) * 8640000
    wdo = np.nanmean([np.nanmean(np.where(sel, divergence(*to_grid(OW[k].real, OW[k].imag), region), np.nan)) for k in range(0, len(OW), 3)]) * 8640000
    print(f'{s:11s} | {fm[2]:8.2f} {fm[3]:10.2f} {fm[0]:8.2f} {fm[1]:+8.1f} {fm[4]:5.2f} |     {fo[2]:6.2f} {fo[3]:6.2f} {fo[0]:7.2f} {fo[1]:+6.1f} {fo[4]:5.2f} | '
          f'{mwe.real:+5.1f}/{mwe.imag:+5.1f} {owe.real:+5.1f}/{owe.imag:+5.1f} | {wdm:+7.1f} {wdo:+7.1f}')

# ---- 1. pressure pattern ----
msl = xr.open_mfdataset(sorted(f for f in glob.glob(f'{oi}/atm_remapped_1m_msl_*.nc') if y0 <= int(f[-12:-8]) <= y1), combine='by_coords', use_cftime=True)['msl']
mt = 'time_counter' if 'time_counter' in msl.dims else 'time'
mm = msl.sel({mt: msl[mt].dt.month.isin(MON)}).mean(mt).load() / 100
ob = xr.open_dataset(f'{NC}/mslp.mon.mean.nc')['mslp']
ob = ob.sel(time=slice('2005', '2014')); ob = ob.where(ob.time.dt.month.isin(MON), drop=True).mean('time')
ob = ob / 100 if float(ob.max()) > 2000 else ob
def asl(p):
    p = p.assign_coords(lon=p.lon % 360).sortby('lon').sortby('lat')
    box = p.sel(lat=slice(-75, -60), lon=slice(170, 298))
    k = np.unravel_index(int(np.nanargmin(box.values)), box.shape)
    sect = p.sel(lat=slice(-75, -60), lon=slice(170, 298)).mean()
    z60 = p.sel(lat=slice(-62, -58)).mean()
    return float(box.values[k]), float(box.lat[k[0]]), float(box.lon[k[1]]) - 360, float(box.values[k] - sect), float(z60)
print(f'\nJJA Amundsen Sea Low (60-75S, 170E-62W): central hPa, lat, lon, relative central pressure (min - sector mean), zonal-mean 60S')
for n, p in (('model', mm), ('NCEP-2 2005-14', ob)):
    a = asl(p); print(f'  {n:16s} {a[0]:7.1f} {a[1]:6.1f} {a[2]:7.1f} {a[3]:+6.1f}   60S mean {a[4]:7.1f}')
