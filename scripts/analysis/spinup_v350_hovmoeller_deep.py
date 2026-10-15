"""Hovmoeller of the global-mean ocean temperature and salinity of the v3.5.0-cc spin-up against PHC3,
and the deep ocean by basin.

Per year: volume-weighted profile of potential temperature and salinity minus PHC3 annual (Steele et al.
2001; in-situ temperature converted to potential with gsw, nearest PHC3 point and level per node), open
ocean only. Per basin and depth band (same basins and bands as deep_warming_where.py): bias of the last
five years against PHC3 and the trend (last five minus years 6-10) per decade.

Usage: spinup_v350_hovmoeller_deep.py <exp root> <y0> <y1> <tag> [nproc]
(module analysis-toolbox/python-04.2026; writes data/clim/hovm_<tag>.npz, data/clim/deep_by_basin_<tag>.txt, plots/<tag>/hovmoeller_temp_salt.png)
"""
import os, sys
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import xarray as xr
import gsw
from multiprocessing import Pool

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
ROOT, Y0, Y1, TAG = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
NP = int(sys.argv[5]) if len(sys.argv) > 5 else 8
BANDS = [(200, 500), (500, 1500), (1500, 2500), (2500, 4000), (4000, 7000)]
md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lat = md['lat'].values; lon = ((md['lon'].values + 180) % 360) - 180
zi = np.abs(md['nz'].values); zm, dz = 0.5 * (zi[1:] + zi[:-1]), np.diff(zi)
vol = md['nod_area'].values[:len(zm)].T * dz[None, :]
open_oc = md['ulevels_nod2D'].values == 1
mid = (lat >= -35) & (lat < 65)
atl = mid & (((lat < 10) & (lon > -70) & (lon < 20)) | ((lat >= 10) & (lon > -98) & (lon < 20) & ~((lat < 18) & (lon < -84))) | ((lat >= 30) & (lon >= 20) & (lon < 45)))
ind = mid & ~atl & (lon >= 20) & (lon < 120) & (lat < 30)
REG = [('Antarctic <60S', lat < -60), ('Southern 60-35S', (lat >= -60) & (lat < -35)), ('Atlantic', atl), ('Indian', ind),
       ('Pacific', mid & ~atl & ~ind), ('Arctic+Nordic >65N', lat >= 65), ('global', np.ones_like(lat, bool))]
ph = xr.open_dataset(f'{P}/input/fesom2/hydrography_recom/phc3.0_annual.nc')
dep = ph['depth'].values[:, None, None]
pt = gsw.pt0_from_t(gsw.SA_from_SP(ph['salt'].values, dep, 0, 0), ph['temp'].values, dep)
iy = np.abs(ph['lat'].values[None, :] - lat[:, None]).argmin(1)
ix = np.abs(((ph['lon'].values[None, :] - lon[:, None] + 180) % 360) - 180).argmin(1)
iz = np.abs(ph['depth'].values[None, :] - zm[:, None]).argmin(1)
PT = pt[iz[None, :], iy[:, None], ix[:, None]]; PS = ph['salt'].values[iz[None, :], iy[:, None], ix[:, None]]


def year(y):
    T = xr.open_dataset(f'{ROOT}/outdata/fesom/temp.fesom.{y}.nc')['temp'].mean('time').values
    S = xr.open_dataset(f'{ROOT}/outdata/fesom/salt.fesom.{y}.nc')['salt'].mean('time').values
    ok = np.isfinite(S) & (S > 1) & open_oc[:, None] & np.isfinite(PT) & np.isfinite(PS)
    V = np.where(ok, vol, 0.0)
    prof = lambda x: (np.where(ok, x, 0) * V).sum(0) / np.maximum(V.sum(0), 1)
    bas = np.full((len(REG), len(BANDS), 2), np.nan)
    for i, (_, reg) in enumerate(REG):
        for j, (a, b) in enumerate(BANDS):
            k = (zm >= a) & (zm < b); w = V[reg][:, k]
            if w.sum() > 0:
                bas[i, j, 0] = (np.where(ok, T - PT, 0)[reg][:, k] * w).sum() / w.sum()
                bas[i, j, 1] = w.sum()
    return prof(T - PT), prof(S - PS), prof(T), bas


if __name__ == '__main__':
    ys = list(range(Y0, Y1 + 1))
    with Pool(NP) as p:
        res = p.map(year, ys)
    dT = np.array([r[0] for r in res]); dS = np.array([r[1] for r in res]); Tm = np.array([r[2] for r in res]); bas = np.array([r[3] for r in res])
    np.savez(f'{REPO}/data/clim/hovm_{TAG}.npz', years=ys, z=zm, dT=dT, dS=dS, T=Tm, basin=bas)
    A, B = slice(5, 10), slice(-5, None); ny = (np.mean(ys[-5:]) - np.mean(ys[5:10])) / 10
    out = [f'{TAG}: potential temperature by basin and depth. bias = {ys[-5]}-{ys[-1]} minus PHC3 annual [K];',
           f'trend = ({ys[-5]}-{ys[-1]} minus {ys[5]}-{ys[9]}) per decade [mK]; first = {ys[0]}-{ys[4]} minus PHC3 [K].', '',
           f"{'basin':20s}" + ''.join(f"{f'{a}-{b} m':>24s}" for a, b in BANDS), f"{'':20s}" + ''.join(f"{'first  bias  trend':>24s}" for _ in BANDS)]
    for i, (name, _) in enumerate(REG):
        line = f'{name:20s}'
        for j in range(len(BANDS)):
            x = bas[:, i, j, 0]
            line += f'{x[:5].mean():+11.2f} {x[B].mean():+5.2f} {1000 * (x[B].mean() - x[A].mean()) / ny:+6.1f}'
        out.append(line)
    out += ['', 'global-mean profile, mK per decade over the same two periods, by level:',
            '  '.join(f'{z:.0f}m:{1000 * (dT[B].mean(0)[k] - dT[A].mean(0)[k]) / ny:+.0f}' for k, z in enumerate(zm) if k % 3 == 0)]
    open(f'{REPO}/data/clim/deep_by_basin_{TAG}.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    os.makedirs(f'{REPO}/plots/{TAG}', exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    for c, (d, lab, vm) in enumerate(((dT, 'potential temperature minus PHC3 (K)', 1.0), (dS, 'salinity minus PHC3', 0.2))):
        for r, (za, zb) in enumerate(((0, 1000), (1000, 6000))):
            k = (zm >= za) & (zm <= zb) & np.isfinite(d[0]); v = vm if r == 0 else vm / 4
            m = ax[r, c].pcolormesh(ys, zm[k], d[:, k].T, cmap='RdBu_r', vmin=-v, vmax=v, shading='nearest')
            ax[r, c].set_ylim(zb, za); ax[r, c].set_ylabel('depth (m)'); fig.colorbar(m, ax=ax[r, c])
            ax[r, c].set_title(f'{lab}, {za}-{zb} m' if r == 0 else f'{za}-{zb} m (colour range / 4)', fontsize=10)
    for a in ax[1]: a.set_xlabel('model year')
    fig.suptitle(f'{TAG}: global-mean ocean drift against PHC3 annual'); fig.tight_layout()
    fig.savefig(f'{REPO}/plots/{TAG}/hovmoeller_temp_salt.png', dpi=130)
