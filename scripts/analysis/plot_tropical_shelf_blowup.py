"""The corrected-CORE3 tropical shelf failure, on the native mesh.

WHAT THIS SHOWS.  11X died at OpenIFS step 9802 on a FESOM MPI_ABORT: temperature below
-5 C in the Timor Sea.  It is not one node.  The same disease is already present through
the whole of core3y1, the one-year new-mesh run, so this figure uses core3y1 (12 months
available) against 11W (old mesh, healthy) at the same month.

THE POINT OF THE THIRD COLUMN.  The cold nodes are the SHALLOW ones.  Outlining the
sub-10 C nodes on top of bottom depth is the whole diagnosis in one panel: the failure
follows bathymetry, not geography, and it is confined to shelves.

MESH DISCIPLINE.  core3y1 is the new mesh (220509 nodes), 11W the old one (211567).  Each
row is drawn on its OWN triangulation; the two are not interchangeable and meshguard
asserts the pairing rather than letting numpy broadcast silently.

SIGN OF fh.  As stored, FESOM's fh is positive UPWARD, out of the ocean.  Verified against
the physical anchor that the tropical ocean must gain heat in the annual mean: taken as
written it makes the tropics lose 13.7 W/m2 and the Southern Ocean gain 18.6, which is
backwards.  Panels are therefore labelled "ocean heat LOSS" for positive.
"""
import os, sys, glob
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import meshguard as mg

OUT = ('/work/ab0246/a270092/postprocessing/'
       'investigation_awiesm3_high_lat_cold_bias_global_toa_positive/'
       'tropical_shelf_blowup_core3.png')
LAT0, LAT1 = -35.0, 35.0
MONTH = 11                     # December, where core3y1's disease is fully developed
BLOWUP = (133.03, -9.86)       # the node FESOM aborted on

RUNS = [('core3y1  NEW mesh (corrected CORE3)',
         '/work/bb1469/a270092/runtime/awiesm3-develop/core3y1/outdata/fesom', 1959),
        ('11W  OLD mesh (core3_beta), healthy',
         '/work/bb1469/a270092/runtime/awiesm3-v3.4/11W/outdata/fesom', 1350)]


def load(d, var, year):
    f = f'{d}/{var}.fesom.{year}.nc'
    if not os.path.exists(f):
        return None
    with xr.open_dataset(f, decode_times=False) as x:
        n = var if var in x.data_vars else [v for v in x.data_vars if 'time' not in v][-1]
        a = x[n].values
    return a[MONTH] if a.ndim == 2 else a


def mesh_of(nnodes):
    m = mg.mesh_for(nnodes)
    with xr.open_dataset(f'{m}/mesh.nc', decode_times=False) as d:
        lon = d['lon'].values.astype(float)
        lat = d['lat'].values.astype(float)
        tri = d['triag_nodes'].values.astype(int) - 1
        zlev = d['depth'].values.astype(float)
        dlev = d['depth_lev'].values.astype(int)
    bottom = zlev[np.clip(dlev - 1, 0, len(zlev) - 1)]
    return m, lon, lat, tri, bottom


def tropical_tris(lon, lat, tri):
    """Triangles inside the band, excluding those wrapped across the dateline."""
    tl = lat[tri]
    keep = (tl.max(axis=1) > LAT0) & (tl.min(axis=1) < LAT1)
    tlo = lon[tri]
    keep &= (tlo.max(axis=1) - tlo.min(axis=1)) < 90.0
    return tri[keep]


fig, axes = plt.subplots(2, 3, figsize=(19.5, 8.4))
panel = iter('abcdef')

for r, (title, ddir, year) in enumerate(RUNS):
    sst = load(ddir, 'sst', year)
    fh = load(ddir, 'fh', year)
    if sst is None or fh is None:
        continue
    m, lon, lat, tri, bottom = mesh_of(sst.size)
    mg.check(m, sst.size, title)
    t = tropical_tris(lon, lat, tri)
    verts = np.stack([lon[t], lat[t]], axis=-1)
    cold = sst < 10.0

    for c, (field, cmap, vmin, vmax, lab) in enumerate([
            (sst, 'RdYlBu_r', -2, 32, 'SST  [$^\\circ$C]'),
            (-fh, 'RdBu_r', -300, 300, 'net surface heat flux  [W m$^{-2}$]\n'
                                       'blue = ocean LOSS'),
            (bottom, 'Blues', 0, 4000, 'bottom depth  [m]')]):
        ax = axes[r, c]
        fv = field[t].mean(axis=1)
        pc = PolyCollection(verts, array=fv, cmap=cmap, edgecolors='none',
                            rasterized=True)
        pc.set_clim(vmin, vmax)
        ax.add_collection(pc)
        if c == 2:   # outline the failing nodes on top of bathymetry
            tc = t[cold[t].sum(axis=1) >= 2]
            if len(tc):
                ax.add_collection(PolyCollection(
                    np.stack([lon[tc], lat[tc]], axis=-1), facecolors='none',
                    edgecolors='magenta', linewidths=0.18, rasterized=True))
        ax.plot(*BLOWUP, marker='*', ms=15, mfc='yellow', mec='k', mew=0.8, zorder=6)
        ax.set_xlim(-180, 180); ax.set_ylim(LAT0, LAT1)
        ax.set_facecolor('0.85')
        ax.set_title(f'({next(panel)}) {lab.splitlines()[0]}', fontsize=9.5)
        if c == 0:
            ax.set_ylabel(title, fontsize=9)
        cb = fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.13,
                          fraction=0.055, aspect=38)
        cb.set_label(lab, fontsize=8)
        cb.ax.tick_params(labelsize=7)
        ax.tick_params(labelsize=7)
    n_cold = int((cold & (lat >= -30) & (lat < 30)).sum())
    med = float(np.median(bottom[cold & (lat >= -30) & (lat < 30)])) if n_cold else np.nan
    axes[r, 0].text(0.01, 0.03, f'{n_cold} tropical nodes < 10 $^\\circ$C'
                    + (f', median depth {med:.0f} m' if n_cold else ''),
                    transform=axes[r, 0].transAxes, fontsize=8,
                    bbox=dict(fc='white', alpha=0.85, ec='none'))

fig.suptitle('Corrected CORE3 cools the tropical shelves and applies the wrong surface heat flux there.\n'
             'Top: new mesh, December of its first year. Bottom: old mesh, same month, unchanged physics. '
             'Magenta outlines nodes below 10 $^\\circ$C; star marks the node FESOM aborted on.',
             fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(OUT, dpi=145)
print('wrote', OUT)
