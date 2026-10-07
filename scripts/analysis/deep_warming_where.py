"""Where the deep-ocean warming of PICAL_crunveg_tke_albsn082 sits, and how far each basin is from PHC3.

Per basin and depth band, volume-weighted on the CORE3 mesh:
  bias   mean potential temperature of 2205-2209 minus PHC3 annual (Steele et al. 2001; its in-situ
         temperature converted to potential with gsw, nearest PHC3 grid point and level per node)
  trend  (mean 2205-2209 minus mean 2171-2175) per decade, mK
  share  the band's part of the global heat gain of that depth band (volume x trend)
Basins: Antarctic south of 60S; Southern 60-35S; Atlantic, Indian, Pacific north of 35S to 65N;
Arctic and Nordic Seas north of 65N. Cavities excluded.

Usage: deep_warming_where.py     (reval environment, ~30 GB; writes data/clim/deep_warming_where.txt)
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import xarray as xr
import gsw

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
D = f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082/outdata/fesom'
A, B = (2171, 2175), (2205, 2209)
BANDS = [(200, 500), (500, 1500), (1500, 2500), (2500, 4000), (4000, 7000)]

md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lat = md['lat'].values; lon = ((md['lon'].values + 180) % 360) - 180
zi = np.abs(md['nz'].values); zm, dz = 0.5 * (zi[1:] + zi[:-1]), np.diff(zi)
vol = md['nod_area'].values[:len(zm)].T * dz[None, :]                    # (nod, nz1)
open_oc = md['ulevels_nod2D'].values == 1
mid = (lat >= -35) & (lat < 65)
atl = mid & (((lat < 10) & (lon > -70) & (lon < 20)) | ((lat >= 10) & (lon > -98) & (lon < 20) & ~((lat < 18) & (lon < -84))) | ((lat >= 30) & (lon >= 20) & (lon < 45)))
ind = mid & ~atl & (lon >= 20) & (lon < 120) & (lat < 30)
REG = [('Antarctic <60S', lat < -60), ('Southern 60-35S', (lat >= -60) & (lat < -35)), ('Atlantic', atl), ('Indian', ind),
       ('Pacific', mid & ~atl & ~ind), ('Arctic+Nordic >65N', lat >= 65), ('global', np.ones_like(lat, bool))]


def mean(y0, y1):
    acc = 0
    for y in range(y0, y1 + 1):
        acc = acc + xr.open_dataset(f'{D}/temp.fesom.{y}.nc')['temp'].mean('time').values
    return acc / (y1 - y0 + 1)


TA, TB = mean(*A), mean(*B)
wet = np.isfinite(TB) & (TB != 0) & open_oc[:, None]
ph = xr.open_dataset(f'{P}/input/fesom2/hydrography_recom/phc3.0_annual.nc')
pt = gsw.pt0_from_t(gsw.SA_from_SP(ph['salt'].values, ph['depth'].values[:, None, None], 0, 0), ph['temp'].values, ph['depth'].values[:, None, None])
iy = np.abs(ph['lat'].values[None, :] - lat[:, None]).argmin(1)
ix = np.abs(((ph['lon'].values[None, :] - lon[:, None] + 180) % 360) - 180).argmin(1)
iz = np.abs(ph['depth'].values[None, :] - zm[:, None]).argmin(1)
PH = pt[iz[None, :], iy[:, None], ix[:, None]]                            # (nod, nz1)
okp = wet & np.isfinite(PH)
ny = (B[0] + B[1]) / 2 - (A[0] + A[1]) / 2
out = [f'PICAL_crunveg_tke_albsn082: potential temperature by basin and depth. bias = {B[0]}-{B[1]} minus PHC3 annual [K];',
       f'trend = ({B[0]}-{B[1]} minus {A[0]}-{A[1]}) per decade [mK]; share = part of the global heat gain of the band [%].', '',
       f"{'basin':20s}" + ''.join(f"{f'{a}-{b} m':>24s}" for a, b in BANDS), f"{'':20s}" + ''.join(f"{'bias  trend share':>24s}" for _ in BANDS)]
for name, reg in REG:
    line = f'{name:20s}'
    for a, b in BANDS:
        lev = (zm >= a) & (zm < b)
        w = np.where(wet & reg[:, None] & lev[None, :], vol, 0.0); wg = np.where(wet & lev[None, :], vol, 0.0)
        wp = np.where(okp & reg[:, None] & lev[None, :], vol, 0.0)
        d = np.nan_to_num(TB - TA)
        tr = 1e4 * np.sum(d * w) / w.sum() / ny
        line += f"{np.sum(np.nan_to_num(TB - PH) * wp) / wp.sum():+10.2f}{tr:+7.1f}{100 * np.sum(d * w) / np.sum(d * wg):+7.0f}"
    out.append(line)
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/deep_warming_where.txt', 'w').write(txt + '\n')
