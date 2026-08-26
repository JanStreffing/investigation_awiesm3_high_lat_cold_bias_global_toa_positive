"""Where does the coupled LAND cold bias sit -- which season, which latitude band?

WHY.  11P (constant 1990, so roughly period-matched to ERA5) has a land T2m bias of
-1.41 K against an ocean bias of -0.69 K.  The land-minus-ocean gap of ~0.7 K is the
robust part: a forcing-period offset moves both together, so it survives the caveat that
ERA5 is a 1940-2014 reference.

But Siberian JJA in 11P is essentially spot on (10.649 against ERA5's 10.666), so the
land deficit is NOT the boreal-summer problem the campaign was built around.  A land
bias that is present in the annual mean and absent in boreal summer is a COLD-SEASON
signal, and possibly not a high-latitude one at all.  This resolves it by season and
band so the next lever is chosen against a target rather than a hunch.

ERA5 is regridded onto the model grid by linear interpolation.  The model land-sea mask
is used for both, so "land" means the same cells in both fields.

READ THE SEASONS, NOT THE ANNUAL MEAN.  The campaign has twice promoted a lever on an
annual or single-season number that the other seasons contradicted, and the detection
thresholds are per-season for that reason.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
ERA5 = '/work/ab0246/a270092/obs/era5/netcdf/T2M.nc'
Y0, Y1 = 1380, 1389
SEASONS = [('DJF', [12, 1, 2]), ('MAM', [3, 4, 5]),
           ('JJA', [6, 7, 8]), ('SON', [9, 10, 11]), ('ANN', list(range(1, 13)))]
BANDS = [('60-90N', 60, 90), ('30-60N', 30, 60), ('0-30N', 0, 30),
         ('30S-0', -30, 0), ('60-30S', -60, -30)]
ARMS = [('11G', f'{R}/Tuning_test_11G_inppmin50k'),
        ('11N', f'{R}/11N'),
        ('11P', f'{R}/11P')]


def grid():
    f = glob.glob(f'{R}/11P/outdata/oifs/atm_remapped_1m_lsm_*.nc')[0]
    with xr.open_dataset(f, decode_times=False) as d:
        lsm = np.asarray(d['lsm'].values, dtype=float)
        if lsm.ndim == 3:
            lsm = lsm[0]
        return lsm, d['lat'].values, d['lon'].values


def era5_on(mlat, mlon):
    """12-month climatology on the model grid.

    The file is a monthly TIME SERIES -- 300 steps, 1990-2014, "hours since
    1900-01-01" -- not a 12-month climatology.  An earlier version tested for a
    length-12 axis, fell through when that failed, and averaged all 300 months into one
    annual field, then differenced every model season against it.  That produced -22 K
    in DJF and +15 K in JJA at 60-90N with the sign mirrored in the SH: a broken
    reference, not a model bias.  Decode the calendar and group by month.
    """
    with xr.open_dataset(ERA5, decode_times=True) as d:
        da = d['t2m']
        clim = da.groupby('time.month').mean('time')
        la = [c for c in clim.dims if 'lat' in c][0]
        lo = [c for c in clim.dims if 'lon' in c][0]
        out = clim.interp({la: ('lat', mlat), lo: ('lon', mlon % 360)}).values
    out = np.asarray(out, dtype=float)
    out = out - 273.15 if np.nanmean(out) > 100 else out
    # SANITY CHECK: the NH extratropical cycle must peak in July (index 6).
    w = np.cos(np.deg2rad(mlat)); nh = mlat > 30
    cyc = [float(np.average(out[i][nh].mean(axis=1), weights=w[nh])) for i in range(12)]
    pk = int(np.argmax(cyc))
    if pk not in (5, 6, 7):
        raise SystemExit(f'ERA5 month axis wrong: NH cycle peaks at index {pk}, '
                         f'expected Jun-Aug. cycle={np.round(cyc,1)}')
    print(f'  ERA5 check: NH>30N cycle peaks at month index {pk} (Jul=6). OK')
    return out


def model_months(root):
    v = []
    for y in range(Y0, Y1 + 1):
        f = f'{root}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc'
        if glob.glob(f):
            with xr.open_dataset(f, decode_times=False) as d:
                v.append(np.asarray(d['2t'].values, dtype=float))
    return np.mean(v, axis=0) - 273.15 if v else None


def main():
    print(__doc__)
    print('=' * 92)
    lsm, mlat, mlon = grid()
    era = era5_on(mlat, mlon)
    w = np.broadcast_to(np.cos(np.deg2rad(mlat))[:, None], lsm.shape)

    def mean(a2d, sel):
        ww = w * sel * np.isfinite(a2d)
        return float(np.nansum(a2d * ww) / ww.sum()) if ww.sum() else np.nan

    for tag, root in ARMS:
        m = model_months(root)
        if m is None:
            print(f'{tag}: no output'); continue
        print(f'\n### {tag}   LAND T2m bias vs ERA5 [K], {Y0}-{Y1}\n')
        hdr = f'  {"band":8s}' + ''.join(f'{s:>8s}' for s, _ in SEASONS)
        print(hdr)
        for nm, lo_, hi in BANDS:
            band = np.zeros_like(lsm, bool)
            band[(mlat >= lo_) & (mlat < hi), :] = True
            sel = band & (lsm >= 0.5)
            if sel.sum() == 0:
                continue
            row = f'  {nm:8s}'
            for _, months in SEASONS:
                idx = [x - 1 for x in months]
                row += f'{mean(m[idx].mean(axis=0), sel) - mean(era[idx].mean(axis=0), sel):8.2f}'
            print(row)
        allland = (lsm >= 0.5)
        row = f'  {"ALL LAND":8s}'
        for _, months in SEASONS:
            idx = [x - 1 for x in months]
            row += f'{mean(m[idx].mean(axis=0), allland) - mean(era[idx].mean(axis=0), allland):8.2f}'
        print(row)
        row = f'  {"ocean":8s}'
        for _, months in SEASONS:
            idx = [x - 1 for x in months]
            row += f'{mean(m[idx].mean(axis=0), lsm < 0.5) - mean(era[idx].mean(axis=0), lsm < 0.5):8.2f}'
        print(row)


if __name__ == '__main__':
    main()
