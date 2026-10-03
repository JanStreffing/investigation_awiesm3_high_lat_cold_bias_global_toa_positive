"""Annual climate series for a coupled run: T2m by region, OpenIFS high-veg cover, and
FESOM's own hemispheric sea-ice integrals.

WHY THIS EXISTS.  To put the albedo legs (2100+, merged LPJ-GUESS canopy work) on one
axis with the same line before 2100 (levante, mirrored), so a step at 2100 can be told
apart from the cold drift the line was already on (-0.214 K/decade global, -0.511 at
60-90N over 2060-2099).  Reads the remapped monthly OIFS files (regular 400x192 grid,
cos-lat weights) and sivoln/sivols/siextentn/siextents (km3 / km2, computed on the
native mesh inside FESOM, so they cannot pick up the wrong mesh).

Usage:  climate_annual_series.py <tag> <exp root> <y0> <y1> <out csv>
"""
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

SIB = (55.0, 75.0, 60.0, 180.0)


def oifs(root, var, yr):
    fn = f'{root}/outdata/oifs/atm_remapped_1m_{var}_{yr}-{yr}.nc'
    if not os.path.exists(fn):
        return None
    with xr.open_dataset(fn) as ds:
        return ds[var].load()


def fesom(root, var, yr):
    fn = f'{root}/outdata/fesom/{var}.fesom.{yr}.nc'
    if not os.path.exists(fn):
        return None
    with xr.open_dataset(fn) as ds:
        return ds[var].values.astype(float)


def wmean(x, w, mask=None):
    w = w if mask is None else w.where(mask, 0.0)
    return float((x * w).sum(('lat', 'lon')) / w.sum(('lat', 'lon')))


def year_row(root, yr, lsm):
    t = oifs(root, '2t', yr)
    if t is None:
        return None
    w = np.cos(np.deg2rad(t.lat)) * xr.ones_like(t.lon)
    lon = t.lon % 360
    land = lsm > 0.5
    sib = land & (t.lat >= SIB[0]) & (t.lat <= SIB[1]) & (lon >= SIB[2]) & (lon <= SIB[3])
    ann = t.mean('time_counter') - 273.15
    djf = t.isel(time_counter=[0, 1, 11]).mean('time_counter') - 273.15
    jja = t.isel(time_counter=[5, 6, 7]).mean('time_counter') - 273.15
    r = {'year': yr,
         't2m_glob': wmean(ann, w),
         't2m_nh': wmean(ann, w, t.lat >= 0),
         't2m_sh': wmean(ann, w, t.lat < 0),
         't2m_6090n': wmean(ann, w, t.lat >= 60),
         't2m_6090s': wmean(ann, w, t.lat <= -60),
         't2m_land': wmean(ann, w, land),
         't2m_ocean': wmean(ann, w, ~land),
         't2m_4590n_land': wmean(ann, w, land & (t.lat >= 45)),
         't2m_sib': wmean(ann, w, sib),
         't2m_sib_djf': wmean(djf, w, sib),
         't2m_sib_jja': wmean(jja, w, sib),
         't2m_6090n_djf': wmean(djf, w, t.lat >= 60),
         't2m_6090n_jja': wmean(jja, w, t.lat >= 60)}
    cvh = oifs(root, 'cvh', yr)
    if cvh is not None:
        r['cvh_sib'] = wmean(cvh.mean('time_counter'), w, sib)
        r['cvh_4590n_land'] = wmean(cvh.mean('time_counter'), w, land & (t.lat >= 45))
    for v in ('tsr', 'ttr'):
        x = oifs(root, v, yr)
        if x is not None:
            r[f'{v}_raw'] = wmean(x.mean('time_counter'), w)
    for v, m in (('sivoln', 3), ('sivoln', 8), ('sivols', 1), ('sivols', 8),
                 ('siextentn', 2), ('siextentn', 8), ('siextents', 1), ('siextents', 8)):
        x = fesom(root, v, yr)
        if x is not None and len(x) == 12:
            r[f'{v}_m{m + 1:02d}'] = x[m]
            r[f'{v}_ann'] = x.mean()
    return r


def main():
    if len(sys.argv) != 6:
        print(__doc__)
        sys.exit(2)
    tag, root, y0, y1, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    lsm = None
    for yr in range(y0, y1 + 1):
        lsm = oifs(root, 'lsm', yr)
        if lsm is not None:
            lsm = lsm.isel(time_counter=0)
            break
    if lsm is None:
        # the mirror carries no lsm: take it from the albedo run, same grid
        with xr.open_dataset('/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4/'
                             'PICAL_ccnice/outdata/oifs/atm_remapped_1m_lsm_2110-2110.nc') as ds:
            lsm = ds['lsm'].isel(time_counter=0).load()
    rows = [r for yr in range(y0, y1 + 1) if (r := year_row(root, yr, lsm)) is not None]
    df = pd.DataFrame(rows)
    df.insert(0, 'exp', tag)
    df.to_csv(out, index=False, float_format='%.6g')
    print(df[['year', 't2m_glob', 't2m_6090n', 't2m_sib', 'sivoln_m04', 'sivols_m09']]
          .to_string(index=False, float_format='%.3f'))


if __name__ == '__main__':
    main()
