"""Southern Ocean surface freshwater budget from FESOM's fw on nodes: cavity melt and open ocean.

fw is the total freshwater flux into the ocean [m/s, positive into the ocean].  Under ice-shelf
cavities (ulevels_nod2D > 1) it is the basal melt; observed Antarctic basal melt is about
1100-1300 Gt/yr (Rignot et al. 2013, Adusumilli et al. 2020).  Over the open ocean 60-78S it is
P-E + runoff + sea-ice growth/melt.  Also the sea-ice part (fw_ice) where the run writes it.
Annual means over Y0-Y1, area-integrated in Gt/yr and as mm/yr per unit area.

Usage:  ROOT=... ARM=PI Y0=2000 Y1=2014 MESH=<mesh diag> python3 scripts/analysis/so_freshwater_budget_nodes.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ['ROOT']; ARM = os.environ['ARM']; MESH = os.environ['MESH']; Y0, Y1 = int(os.environ['Y0']), int(os.environ['Y1'])
with xr.open_dataset(MESH) as m:
    lat = m['lat'].values; ul = m['ulevels_nod2D'].values if 'ulevels_nod2D' in m else np.ones(len(lat), int)
    area = m['nod_area'].values[ul - 1, np.arange(len(lat))]          # area at each node's top wet level (under cavities that is not level 1)
cav = ul > 1; so = (lat >= -78) & (lat <= -60) & ~cav; sh = (lat < -60) & ~cav
GT = 1000.0 * 365.25 * 86400 / 1e12                                   # m/s * m2 -> Gt/yr


def ann(v):
    acc = 0; n = 0
    for y in range(Y0, Y1 + 1):
        f = f'{R}/{ARM}/outdata/fesom/{v}.fesom.{y}.nc'
        if not os.path.exists(f): continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = np.squeeze(d[v].values).astype('f8')
        a[np.abs(a) > 1e10] = np.nan; acc = acc + np.nanmean(a, 0); n += 1
    return (acc / n) if n else None


fw = ann('fw'); fwi = ann('fw_ice'); ro = ann('runoff')
print(f'{ARM} {Y0}-{Y1}: freshwater into the ocean (positive = freshening)')
print(f'  cavity nodes: {cav.sum()} ({area[cav].sum() / 1e12:.2f} M km2); fw under cavities = basal melt {np.nansum(fw[cav] * area[cav]) * GT:8.0f} Gt/yr'
      f'  ({np.nansum(fw[cav] * area[cav]) / area[cav].sum() * 365.25 * 86400 * 1000:.0f} mm/yr)')
for name, k in (('open ocean 60-78S', so), ('open ocean south of 60S', sh)):
    print(f'  {name}: total fw {np.nansum(fw[k] * area[k]) * GT:8.0f} Gt/yr ({np.nansum(fw[k] * area[k]) / area[k].sum() * 365.25 * 86400 * 1000:+.0f} mm/yr)'
          + (f', of which sea ice {np.nansum(fwi[k] * area[k]) * GT:+.0f} Gt/yr' if fwi is not None else '')
          + (f', runoff+calving {np.nansum(ro[k] * area[k]) * GT:+.0f} Gt/yr' if ro is not None else ''))
print('  (sign as written by FESOM; see the source check for the convention)')
