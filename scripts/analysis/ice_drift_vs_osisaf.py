"""Is the model's winter ice drift and divergence the right size? FESOM daily ice velocity against the
OSI SAF sea-ice drift climate data record (OSI-455 v1.0, merged, 24 h displacements, EASE2 75 km LAEA).

Observations: satellite vectors only (status_flag 20, 21, 22, 30); wind-filled or blended vectors
(23-25) are dropped so the reference does not contain a wind model. Speed = |dX,dY| / (t1 - t0).
Model: daily uice/vice (east/north, 0.5 deg regular FESOM output) interpolated bilinearly to the OSI
grid centres and rotated into the grid's x/y axes with pyproj, daily a_ice likewise.
Both: divergence du/dx + dv/dy by centred differences on the same 75 km grid, only where all four
neighbours are valid, so both fields see the same stencil and smoothing.
Common region: cells where OSI has a valid vector on >= 50 % of that month's days across the obs
years; the model is additionally required to have a_ice >= 0.15 that day.
Sectors (SH): Weddell 60W-20E, Indian 20E-90E, Pacific 90E-160E, Ross 160E-130W, Amundsen-Bell. 130W-60W.
Usage: python ice_drift_vs_osisaf.py <model outdata/fesom> <y0> <y1> <label> <sh|nh>
"""
import sys, glob, numpy as np, xarray as xr, pyproj, warnings
warnings.filterwarnings('ignore')
d, y0, y1, label, hem = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
OBS = sorted(glob.glob(f'/albedo/work/user/jstreffi/obs/osisaf_drift/{hem}/*.nc'))
MON = [6, 7, 8] if hem == 'sh' else [12, 1, 2]
GOOD = [20, 21, 22, 30]

g0 = xr.open_dataset(OBS[0])
lat, lon = g0['lat'].values, g0['lon'].values
proj = pyproj.Proj(g0['Lambert_Azimuthal_Equal_Area'].attrs['proj4_string'])
# unit vectors of east and north in grid x/y, numerically
x0, y0g = proj(lon, lat)
xe, ye = proj(lon + 0.01, lat); xs, ys = proj(lon, lat - 0.01)   # east step, and a south step (valid at both poles' ice)
ex, ey = xe - x0, ye - y0g; nx, ny = x0 - xs, y0g - ys              # north = minus south
ne = np.hypot(ex, ey); nn = np.hypot(nx, ny)
ex, ey, nx, ny = ex / ne, ey / ne, nx / nn, ny / nn
DX = 75e3

def divergence(u, v, ok):
    """u, v in m/s on the grid axes; returns 1/s, NaN unless all four neighbours are ok."""
    div = np.full(u.shape, np.nan)
    c = ok[1:-1, 1:-1] & ok[1:-1, 2:] & ok[1:-1, :-2] & ok[2:, 1:-1] & ok[:-2, 1:-1]
    # yc decreases with row index in these files: dv/dy = (v[row-1] - v[row+1]) / (2 DX)
    sx = 1.0 if float(g0.xc[1] - g0.xc[0]) > 0 else -1.0
    dudx = sx * (u[1:-1, 2:] - u[1:-1, :-2]) / (2 * DX)
    sgn = 1.0 if float(g0.yc[1] - g0.yc[0]) < 0 else -1.0
    dvdy = sgn * (v[:-2, 1:-1] - v[2:, 1:-1]) / (2 * DX)
    div[1:-1, 1:-1] = np.where(c, dudx + dvdy, np.nan)
    return div

# ---------------- observations ----------------
obs = {m: {'spd': [], 'div': [], 'ok': []} for m in MON}
for f in OBS:
    with xr.open_dataset(f) as o:
        m = int(o.time.dt.month[0])
        if m not in MON: continue
        fl = o['status_flag'].values[0]; ok = np.isin(fl, GOOD)
        dt = (o['t1'].values[0] - o['t0'].values[0]).astype('timedelta64[s]').astype(float)
        dt = np.where(np.isfinite(dt) & (dt > 0), dt, 86400.0)
        u = o['dX'].values[0] * 1e3 / dt; v = o['dY'].values[0] * 1e3 / dt
        ok &= np.isfinite(u) & np.isfinite(v)
        obs[m]['spd'].append(np.where(ok, np.hypot(u, v), np.nan))
        obs[m]['div'].append(divergence(np.where(ok, u, 0), np.where(ok, v, 0), ok))
        obs[m]['ok'].append(ok)
region = {m: np.mean(obs[m]['ok'], 0) >= 0.5 for m in MON}
nyrs = len({f.split('-')[-1][:4] for f in OBS})

# ---------------- model ----------------
def interp_regular(field, glat, glon):
    """bilinear from a regular lat/lon grid (lat ascending, lon 0..360 or -180..180) to points"""
    la = glat; lo = np.mod(glon, 360)
    order = np.argsort(lo); lo = lo[order]; field = field[:, order]
    lo_ext = np.concatenate([lo[-1:] - 360, lo, lo[:1] + 360])
    f_ext = np.concatenate([field[:, -1:], field, field[:, :1]], 1)
    if la[0] > la[-1]: la = la[::-1]; f_ext = f_ext[::-1]
    pl, pn = lat.ravel(), np.mod(lon.ravel(), 360)
    i = np.clip(np.searchsorted(la, pl) - 1, 0, len(la) - 2); j = np.clip(np.searchsorted(lo_ext, pn) - 1, 0, len(lo_ext) - 2)
    wy = (pl - la[i]) / (la[i + 1] - la[i]); wx = (pn - lo_ext[j]) / (lo_ext[j + 1] - lo_ext[j])
    out = (f_ext[i, j] * (1 - wx) * (1 - wy) + f_ext[i, j + 1] * wx * (1 - wy) +
           f_ext[i + 1, j] * (1 - wx) * wy + f_ext[i + 1, j + 1] * wx * wy)
    return out.reshape(lat.shape)

mod = {m: {'spd': [], 'div': []} for m in MON}
for y in range(y0, y1 + 1):
    U = xr.open_dataset(f'{d}/uice.fesom.gr.{y}.nc', decode_times=False)
    V = xr.open_dataset(f'{d}/vice.fesom.gr.{y}.nc', decode_times=False)
    A = xr.open_dataset(f'{d}/a_ice.fesom.gr.{y}.nc', decode_times=False)
    glat, glon = U['lat'].values, U['lon'].values
    nt = U.sizes['time']; doy = np.arange(nt)
    month = np.searchsorted(np.cumsum([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]), doy, side='right') + 1
    for t in np.where(np.isin(month, MON))[0]:
        m = int(month[t])
        ue = np.nan_to_num(U['uice'].values[t]); vn = np.nan_to_num(V['vice'].values[t])
        ai = np.nan_to_num(A['a_ice'].values[t])
        ue_p, vn_p, ai_p = (interp_regular(x, glat, glon) for x in (ue, vn, ai))
        ug = ue_p * ex + vn_p * nx; vg = ue_p * ey + vn_p * ny
        ok = region[m] & (ai_p >= 0.15)
        mod[m]['spd'].append(np.where(ok, np.hypot(ug, vg), np.nan))
        mod[m]['div'].append(divergence(np.where(ok, ug, 0), np.where(ok, vg, 0), ok))

SECT = {'all': np.ones_like(lon, bool)}
if hem == 'sh':
    L = np.mod(lon + 180, 360) - 180
    SECT.update({'Weddell': (L >= -60) & (L < 20), 'Indian': (L >= 20) & (L < 90), 'Pacific': (L >= 90) & (L < 160),
                 'Ross': (L >= 160) | (L < -130), 'Amund-Bell': (L >= -130) & (L < -60)})
print(f'{label} {y0}-{y1} vs OSI-455 ({nyrs} winters, satellite vectors only), {hem.upper()} months {MON}')
print('speed cm/s; divergence in %/day (1e-2 per day = fractional area change per day); |div| = mean absolute')
print(f'{"sector":11s} {"month":>5s} {"spd obs":>8s} {"spd mod":>8s} {"ratio":>6s} {"div obs":>8s} {"div mod":>8s} {"|div| obs":>9s} {"|div| mod":>9s} {"ratio":>6s}')
for s, fn in SECT.items():
    sel = fn
    for m in MON + ['win']:
        ms = MON if m == 'win' else [m]
        def agg(store, key, f):
            arr = np.concatenate([np.stack(store[k][key]) for k in ms])
            mask = np.concatenate([np.broadcast_to(region[k] & sel, np.stack(store[k][key]).shape) for k in ms])
            return f(np.where(mask, arr, np.nan))
        so, sm = agg(obs, 'spd', np.nanmean) * 100, agg(mod, 'spd', np.nanmean) * 100
        do, dm = agg(obs, 'div', np.nanmean) * 86400 * 100, agg(mod, 'div', np.nanmean) * 86400 * 100
        ao, am = (agg(obs, 'div', lambda x: np.nanmean(np.abs(x))) * 86400 * 100,
                  agg(mod, 'div', lambda x: np.nanmean(np.abs(x))) * 86400 * 100)
        if m != 'win' and s != 'all': continue
        print(f'{s:11s} {str(m):>5s} {so:8.2f} {sm:8.2f} {sm/so:6.2f} {do:+8.3f} {dm:+8.3f} {ao:9.3f} {am:9.3f} {am/ao:6.2f}')
