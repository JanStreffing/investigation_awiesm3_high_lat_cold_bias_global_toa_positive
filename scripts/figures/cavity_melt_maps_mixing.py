"""Ice-shelf basal melt under cavities, per vertical-mixing scheme, on a polar projection.

fw (freshwater volume flux, m s-1, positive OUT of the ocean) at cavity nodes is minus the
basal melt, so melt = -fw: positive where the shelf loses mass, NEGATIVE where marine ice
accretes.  That sign change is real physics under the cold cavities, so every colour scale
here is diverging and symmetric about zero rather than a one-sided magnitude ramp.

Cavity nodes are ulevels_nod2D > 1: the ocean top there is the shelf base, not the sea
surface.  Their area must be taken at each node's own first wet level, because nod_area is
identically zero in the surface layer for exactly these nodes.

Needs the conda env with a working cartopy; the spack mambaforge one cannot load shapely:
  /home/a/a270092/.conda/envs/esm-tools_auto_tripyview/bin/python3

Usage:  Y0=2030 Y1=2039 <that python> scripts/figures/cavity_melt_maps_mixing.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings, textwrap
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import matplotlib.path as mpath
from matplotlib.colors import SymLogNorm
import cartopy.crs as ccrs, cartopy.feature as cfeature
warnings.filterwarnings('ignore')

REPO = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive'
R    = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESH = '/work/ab0246/a270092/input/fesom2/core3/fesom.mesh.diag.nc'
Y0, Y1 = int(os.environ.get('Y0', 2030)), int(os.environ.get('Y1', 2039))
RHOW, SEC = 1000.0, 365.25*86400.0
SURF, INK, MUTED = '#fcfcfb', '#0b0b0b', '#52514e'
LAT_N = -60.0

ARMS = [('PICAL_ccnice', 'KPP  (control)'), ('PICAL_cvKPP', 'cvmix_KPP'),
        ('PICAL_cvTKE', 'cvmix_TKE'), ('PICAL_cvTKEIDEMIX', 'cvmix_TKE + cvmix_IDEMIX')]

with xr.open_dataset(MESH) as m:
    lon = np.asarray(m['lon'].values, float)          # already degrees
    lat = np.asarray(m['lat'].values, float)
    ule = np.asarray(m['ulevels_nod2D'].values, int)
    tri = np.asarray(m['face_nodes'].values, int).T - 1   # (3, elem) one-based
    _na = np.asarray(m['nod_area'].values, float)
    narea = _na[ule-1, np.arange(_na.shape[1])]       # area at each node's first wet level
cav = ule > 1

# Project the nodes ourselves and triangulate in projection coordinates.  Handing cartopy a
# Triangulation plus transform=PlateCarree() silently draws nothing, so the map came out
# blank; in projection space tripcolor is plain matplotlib and behaves.
PROJ = ccrs.SouthPolarStereo(central_longitude=0)
_xy = PROJ.transform_points(ccrs.PlateCarree(), lon, lat)
px, py = _xy[:, 0], _xy[:, 1]
# only triangles that lie wholly inside a cavity; longitude wrap is harmless on a polar map
tri_ok = np.all((cav & (lat < LAT_N))[tri], axis=1)
t = Triangulation(px, py, tri[tri_ok])
print(f'cavity nodes {cav.sum()}, drawn triangles {tri_ok.sum()}, cavity area {narea[cav].sum():.3e} m2')

def melt(arm):
    acc, n = None, 0
    for y in range(Y0, Y1+1):
        p = f'{R}/{arm}/outdata/fesom/fw.fesom.{y}.nc'
        if not os.path.exists(p): continue
        with xr.open_dataset(p, decode_times=False) as d:
            a = np.asarray(d['fw'].values, float)
        a = a.mean(0) if a.ndim > 1 else a
        acc = a if acc is None else acc + a; n += 1
    return None if n == 0 else -(acc/n)*SEC           # m/yr water equivalent, melt positive

def integral(f):
    return float(np.nansum(f[cav]*narea[cav])*RHOW/1e12)   # Gt/yr

fields = [melt(a) for a, _ in ARMS]
totals = [integral(f) if f is not None else float('nan') for f in fields]
for (arm, _), tot in zip(ARMS, totals):
    print(f'{arm:<20} total basal melt = {tot:8.1f} Gt/yr')

# Melt spans four decades (median 0.4, peak 40 m/yr) with ~5 % refreezing, so a linear
# scale would render almost every cavity white.  A symmetric-log norm keeps zero at the
# centre of the diverging map and still resolves both the refreezing and the grounding-line
# hotspots.
dif  = [f - fields[0] for f in fields[1:]]
vabs = 40.0
vdif = 20.0
nrm_abs = SymLogNorm(linthresh=0.1, linscale=0.6, vmin=-vabs, vmax=vabs, base=10)
nrm_dif = SymLogNorm(linthresh=0.05, linscale=0.6, vmin=-vdif, vmax=vdif, base=10)
TK_ABS = [-10, -1, -0.1, 0, 0.1, 1, 10]
TK_DIF = [-10, -1, -0.1, 0, 0.1, 1, 10]

proj = PROJ
circle = mpath.Path(np.column_stack([0.5 + 0.5*np.cos(np.linspace(0, 2*np.pi, 200)),
                                     0.5 + 0.5*np.sin(np.linspace(0, 2*np.pi, 200))]))

def polar_ax(ax):
    ax.set_extent([-180, 180, -90, LAT_N], ccrs.PlateCarree())
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.add_feature(cfeature.LAND, facecolor='#e9e7e1', edgecolor='none', zorder=1)
    ax.coastlines(resolution='50m', linewidth=0.45, color='#8d8b85', zorder=5)
    ax.gridlines(linewidth=0.35, color='#c6c4be', alpha=0.9, zorder=4,
                 xlocs=range(-180, 181, 60), ylocs=[-80, -70, -60])
    ax.spines['geo'].set_edgecolor('#b8b6b0'); ax.spines['geo'].set_linewidth(0.6)

fig = plt.figure(figsize=(16.4, 9.6))
fig.patch.set_facecolor(SURF)
gs = fig.add_gridspec(2, 4, left=0.02, right=0.98, top=0.865, bottom=0.115,
                      wspace=0.06, hspace=0.20)

for j, ((arm, lab), f, tot) in enumerate(zip(ARMS, fields, totals)):
    ax = fig.add_subplot(gs[0, j], projection=proj); polar_ax(ax)
    if f is not None:
        im0 = ax.tripcolor(t, f, cmap='RdBu_r', norm=nrm_abs,
                           shading='gouraud', zorder=3, rasterized=True)
    ax.set_title(f'{lab}\n{tot:.0f} Gt/yr', fontsize=11, color=INK, pad=6)

for j, (d, (arm, lab)) in enumerate(zip(dif, ARMS[1:]), start=1):
    ax = fig.add_subplot(gs[1, j], projection=proj); polar_ax(ax)
    im1 = ax.tripcolor(t, d, cmap='BrBG_r', norm=nrm_dif,
                       shading='gouraud', zorder=3, rasterized=True)
    ax.set_title(f'{lab} - control\n{totals[j]-totals[0]:+.0f} Gt/yr', fontsize=11, color=INK, pad=6)

cax0 = fig.add_axes([0.036, 0.360, 0.185, 0.015])
cb0 = fig.colorbar(im0, cax=cax0, orientation='horizontal', ticks=TK_ABS)
cb0.ax.set_xticklabels([f'{v:g}' for v in TK_ABS])
cb0.set_label('basal melt rate  [m/yr w.e.]\n+ melting  /  - refreezing', fontsize=8.5, color=MUTED)
cax1 = fig.add_axes([0.036, 0.235, 0.185, 0.015])
cb1 = fig.colorbar(im1, cax=cax1, orientation='horizontal', ticks=TK_DIF)
cb1.ax.set_xticklabels([f'{v:g}' for v in TK_DIF])
cb1.set_label('change in basal melt rate  [m/yr w.e.]', fontsize=8.5, color=MUTED)
for cb in (cb0, cb1):
    cb.ax.tick_params(labelsize=8, colors=MUTED); cb.outline.set_linewidth(0.4)

fig.text(0.02, 0.970, f'Ice-shelf basal melt under cavities, {Y0}-{Y1}',
         fontsize=15, color=INK, ha='left', va='center')
fig.text(0.02, 0.934, 'absolute rates per mixing scheme (top), and the change each scheme makes '
         'against the KPP control (bottom)', fontsize=10, color=MUTED, ha='left', va='center')
note = ('Symmetric-log scales, zero at the centre of the diverging map, linear within +/-0.1. '
        'Melt = -fw at cavity nodes (ulevels_nod2D > 1), m/yr water equivalent; totals integrate over each node\'s '
        'area at its own first wet level, x rho_water; scales run to +/-40 (top) and +/-20 (bottom) m/yr. '
        'For scale: Rignot 2013 ~1500 Gt/yr, Adusumilli 2020 ~1100 Gt/yr.')
fig.text(0.02, 0.016, '\n'.join(textwrap.wrap(note, 178)), fontsize=8, color=MUTED, va='bottom')
p = f'{REPO}/plots/cavity_melt_maps_mixing.png'
fig.savefig(p, dpi=170, facecolor=SURF); plt.close(fig); print('wrote', p)
