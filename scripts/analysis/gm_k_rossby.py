"""GM coefficient from the local Rossby radius instead of the mesh spacing.

Rd = min(c1/max(|f|,1e-6), 200 km), c1 = max(0.1, sum(N dz)/pi), as FESOM computes it in
oce_fer_gm.F90:init_Redi_GM, here from gm2500's (PICAL_crunveg_ob1200) monthly N2, averaged over
2121-2129 (annual mean of the monthly Rd, and March / September separately).
Candidates: K = K_max * min(1, (Rd/R_ref)^p), p = 1 and 2, with R_ref chosen so that the subtropics
saturate. Compared region by region with CORE2's AWI-CM3 v3.3 field (interpolated to CORE3 nodes,
inverse-distance, 4 neighbours) and with CORE3 production; maps of Rd and of the candidates.
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.tri as mtri, cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm
from netCDF4 import Dataset
from scipy.spatial import cKDTree
exec(open('scripts/figures/gm_k_maps.py').read().split("C2 = mesh(")[0])          # mesh(), k_core2(), k_core3()
C2 = mesh('/albedo/pool/fesom2/core2/fesom.mesh.diag.nc')
DIAG3 = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/input/fesom2/core3/fesom.mesh.diag.nc'
C3 = mesh(DIAG3)
A3, lo3, la3, tri3 = C3
with Dataset(DIAG3) as nc: zb = np.abs(np.asarray(nc.variables['nz'][:], 'f8'))
dz = np.diff(zb)
D = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4/PICAL_crunveg_ob1200/outdata/fesom'
f = np.maximum(np.abs(2 * 7.292e-5 * np.sin(np.deg2rad(la3))), 1e-6)
rd_m = np.zeros((12, len(la3))); n = 0
for y in range(2121, 2130):
    with Dataset(f'{D}/N2.fesom.{y}.nc') as nc:
        for m in range(12):
            n2 = np.asarray(nc.variables['N2'][m], 'f8'); n2[~np.isfinite(n2) | (np.abs(n2) > 1)] = 0.0
            N = np.sqrt(np.maximum(n2, 0.0))
            c1 = np.maximum(0.1, (0.5 * (N[:, :-1] + N[:, 1:]) * dz).sum(1) / np.pi)
            rd_m[m] += np.minimum(c1 / f, 2e5)
    n += 1
rd_m /= n
Rd = rd_m.mean(0) / 1e3                                                   # km, annual mean
np.save('data/clim/rossby_radius_gm2500_2121-2129_km.npy', Rd)

def xyz(lo, la):
    lo, la = np.deg2rad(lo), np.deg2rad(la); return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
d, i = cKDTree(xyz(C2[1], C2[2])).query(xyz(lo3, la3), k=4); w = 1 / np.maximum(d, 1e-9)
k_c2 = (k_core2(C2[0])[i] * w).sum(1) / w.sum(1)

OPTS = {'CORE2 v3.3 (on CORE3)': k_c2, 'CORE3 production': k_core3(A3, la3),
        'Rd, p=1, Rref 30': 2500 * np.minimum(1, Rd / 30), 'Rd, p=2, Rref 30': 2500 * np.minimum(1, (Rd / 30) ** 2),
        'Rd, p=1, Rref 20': 2500 * np.minimum(1, Rd / 20), 'Rd, p=2, Rref 20': 2500 * np.minimum(1, (Rd / 20) ** 2)}
R = {'Labrador': (la3 >= 56) & (la3 < 62) & (lo3 >= -60) & (lo3 < -50),
     'SPNA 45-66N': (la3 >= 45) & (la3 < 66) & (lo3 >= -60) & (lo3 <= -10),
     'Nordic/Arctic >66N': la3 >= 66, 'N Pacific 45-66N': (la3 >= 45) & (la3 < 66) & ((lo3 < -120) | (lo3 > 140)),
     'tropics 10S-10N': np.abs(la3) < 10, 'subtropics 20-40': (np.abs(la3) >= 20) & (np.abs(la3) < 40),
     'SO 40-55S': (la3 < -40) & (la3 >= -55), 'SO 55-65S': (la3 < -55) & (la3 >= -65), 'Antarctic <65S': la3 < -65}
print(f'{"region":20s}{"Rd km":>7s}{"Mar":>6s}{"Sep":>6s}' + ''.join(f'{k:>23s}' for k in OPTS))
for r, msk in R.items():
    wa = A3[msk]
    print(f'{r:20s}{np.average(Rd[msk], weights=wa):7.1f}{np.average(rd_m[2][msk], weights=wa)/1e3:6.1f}'
          f'{np.average(rd_m[8][msk], weights=wa)/1e3:6.1f}' + ''.join(f'{np.average(v[msk], weights=wa):23.0f}' for v in OPTS.values()))
print(f'{"RMS vs CORE2":20s}{"":19s}' + ''.join(f'{np.sqrt(np.average((v - k_c2)**2, weights=A3)):23.0f}' for v in OPTS.values()))

panels = [('Rossby radius Rd (annual mean, gm2500 2121-2129), km', Rd, np.array([0, 2, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200]), 'Purples'),
          ('CORE2, AWI-CM3 v3.3, interpolated to CORE3', k_c2, None, 'Blues'),
          ('Rd-based, K = 2500 min(1, Rd/20 km)', OPTS['Rd, p=1, Rref 20'], None, 'Blues'),
          ('Rd-based, K = 2500 min(1, (Rd/20 km)^2)', OPTS['Rd, p=2, Rref 20'], None, 'Blues')]
fig = plt.figure(figsize=(14, 9.5))
span = lo3[tri3].max(1) - lo3[tri3].min(1)
for j, (title, v, lev, cm) in enumerate(panels):
    ax = fig.add_subplot(2, 2, j + 1, projection=ccrs.Robinson())
    xy = ax.projection.transform_points(ccrs.PlateCarree(), lo3, la3)
    T = mtri.Triangulation(xy[:, 0], xy[:, 1], tri3[span < 180])
    lev = np.arange(0, 3001, 250) if lev is None else lev
    cmap = plt.get_cmap(cm, len(lev) - 1); norm = BoundaryNorm(lev, cmap.N)
    mm = ax.tripcolor(T, np.clip(v, lev[0], lev[-1] - 1e-6), cmap=cmap, norm=norm, shading='gouraud', rasterized=True)
    ax.coastlines(lw=0.3, color='#555555'); ax.set_global(); ax.set_title(title, fontsize=9.5, loc='left')
    cb = fig.colorbar(mm, ax=ax, orientation='horizontal', fraction=0.05, pad=0.03, ticks=lev if j == 0 else lev[::2])
    cb.ax.tick_params(labelsize=7); cb.set_label('km' if j == 0 else 'm$^2$/s, top level', fontsize=8)
fig.tight_layout(); fig.savefig('plots/gm_k_rossby_maps.png', dpi=120, bbox_inches='tight')
print('written')
