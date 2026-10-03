"""Maps of the GM coefficient at the top level (horizontal factor x K_GM_max, before the vertical
scaling and slope tapering), computed from the mesh with FESOM's own formula
(oce_fer_gm.F90:init_Redi_GM, plus the K_GM_max_NH option of feat/gm-hemispheric 636b26aa):
  CORE2, AWI-CM3 v3.3   K_GM_max 3000, resolution scaling order 1, ramp off < 30 km, full > 40 km
  CORE3, production     K_GM_max 2500, resolution scaling order 2, no ramp
  CORE3, gmhemi         as production, K_GM_max 2500 in the SH, 1000 in the NH, linear over 5S-5N
reso = 2 sqrt(A/pi); order 1: min(1, reso/100 km); order 2: min(1, sqrt(2 A / (100 km)^2)).
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.tri as mtri, cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm
from netCDF4 import Dataset

def mesh(f):
    with Dataset(f) as nc:
        A = np.asarray(nc.variables['nod_area'][0], 'f8')
        if 'lon' in nc.variables:      # newer diag layout
            lon = np.asarray(nc.variables['lon'][:], 'f8'); lat = np.asarray(nc.variables['lat'][:], 'f8')
            tri = np.asarray(nc.variables['face_nodes'][:], 'i8').T
        else:                          # older layout: nodes(2, n), elem(3, n)
            xy = np.asarray(nc.variables['nodes'][:], 'f8'); lon, lat = xy[0], xy[1]
            tri = np.asarray(nc.variables['elem'][:], 'i8').T
    if np.abs(lon).max() < 7: lon, lat = np.rad2deg(lon), np.rad2deg(lat)
    lon = ((lon + 180) % 360) - 180
    if tri.min() == 1: tri -= 1
    return A, lon, lat, tri

def k_core2(A):
    reso = 2 * np.sqrt(A / np.pi) / 1e3
    return 3000 * np.minimum(1, reso / 100) * np.clip((reso - 30) / 10, 0, 1)

def k_core3(A, lat, nh=None):
    s = np.minimum(1, np.sqrt(2 * A / 1e10))
    kmax = np.full_like(A, 2500.0)
    if nh is not None:
        w = np.clip((lat + 5) / 10, 0, 1); kmax = (1 - w) * 2500 + w * nh
    return s * kmax

C2 = mesh('/albedo/pool/fesom2/core2/fesom.mesh.diag.nc')
C3 = mesh('/albedo/work/projects/p_awiesm3_cmip7/jstreffi/input/fesom2/core3/fesom.mesh.diag.nc')
cases = [('CORE2, AWI-CM3 v3.3 (K_GM_max 3000, linear resolution scaling, off below 30-40 km)', C2, k_core2(C2[0])),
         ('CORE3, production (K_GM_max 2500, resolution scaling, no ramp)', C3, k_core3(C3[0], C3[2])),
         ('CORE3, gmhemi (2500 in the SH, 1000 in the NH, linear over 5S-5N)', C3, k_core3(C3[0], C3[2], nh=1000.0))]
levels = np.arange(0, 3001, 250)
cmap = plt.get_cmap('Blues', len(levels) - 1)
norm = BoundaryNorm(levels, cmap.N)
fig = plt.figure(figsize=(11, 14))
for i, (title, (A, lon, lat, tri), k) in enumerate(cases):
    ax = fig.add_subplot(3, 1, i + 1, projection=ccrs.Robinson())
    span = lon[tri].max(1) - lon[tri].min(1)                       # drop dateline-wrapping triangles
    xy = ax.projection.transform_points(ccrs.PlateCarree(), lon, lat)  # project first: cartopy cannot
    T = mtri.Triangulation(xy[:, 0], xy[:, 1], tri[span < 180])        # transform a tripcolor
    m = ax.tripcolor(T, k, cmap=cmap, norm=norm, shading='gouraud', rasterized=True)
    ax.coastlines(lw=0.4, color='#555555'); ax.set_global()
    lab = (lat >= 56) & (lat < 62) & (lon >= -60) & (lon < -50); so = (lat < -55) & (lat >= -65)
    ax.set_title(f'{title}\nLabrador {np.average(k[lab], weights=A[lab]):.0f}, SO 55-65S {np.average(k[so], weights=A[so]):.0f} m$^2$/s',
                 fontsize=10, loc='left')
cb = fig.colorbar(m, ax=fig.axes, orientation='horizontal', fraction=0.03, pad=0.03, ticks=levels[::2])
cb.set_label('GM coefficient at the top level, m$^2$/s (before depth scaling and slope tapering)')
fig.savefig('plots/gm_k_maps_core2_core3.png', dpi=130, bbox_inches='tight')
print('written')
