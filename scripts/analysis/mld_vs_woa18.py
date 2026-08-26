"""Mixed layer depth against WOA18: is the model pumping surface heat down too efficiently?

THE HYPOTHESIS.  The model absorbs +0.8 W/m2 at TOA while sitting 1.1 K too cold over
0-100 m and 0.15 K too warm at 700-2000 m against the PHC3 it was initialised from.  Those
are not two faults, they are one: heat is being removed from the surface faster than it
should be, so the surface never warms, never radiates the excess back as OLR, and the
imbalance persists while the interior fills.  Burying it deeper would make the model LOOK
equilibrated for longer while making the cold bias worse.  The fix runs the other way, keep
the heat where it can radiate.

THE DIRECT TEST.  If excess downward transport is the mechanism, the mixed layer should be
too DEEP.  FESOM writes MLD1 and MLD2; WOA18 carries an objectively analysed mixed layer
thickness climatology.  Comparing them says whether the mechanism is real before any mixing
parameter is touched.

WHAT THIS COMPARISON CAN AND CANNOT SAY.  The WOA18 file on disk is FALL ONLY, October to
December, 1981-2010, so it is NH autumn and SH spring.  Mixed layers are deepest in winter
and that is the season that sets how much heat gets subducted, so this is the wrong season
for the strongest statement and the right one only for an indication.  A winter file, or an
Argo MLD climatology with real Southern Ocean coverage, would be needed to settle it.
WOA18 is also 1981-2010 while the model is a 1850 control; mixed layer depth is not strongly
forced on that timescale but the offset is not zero either.

METHOD.  Band means are computed on each grid in its own geometry and compared as band
means, so no regridding is involved and no interpolation error is introduced.  Model months
10-12 only, to match.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESH = '/work/ab0246/a270092/input/fesom2/core3_beta'
WOA = '/work/ab0246/a270092/obs/WOA/woa18_decav81B0_M0216_01.nc'
ARM = os.environ.get('ARM', '11Q')
YEARS = list(range(1382, 1390))
FALL = [9, 10, 11]            # 0-based Oct, Nov, Dec

BANDS = [('90-60S  Antarctic/SO', -90, -60), ('60-45S  subantarctic', -60, -45),
         ('45-30S', -45, -30), ('30S-30N  tropics', -30, 30), ('30-45N', 30, 45),
         ('45-60N  subpolar NA', 45, 60), ('60-90N  Arctic/Nordic', 60, 90)]

root = f'{R092}/{ARM}' if os.path.isdir(f'{R092}/{ARM}') else \
       [p for p in glob.glob(f'{R092}/*{ARM}*') if os.path.isdir(p)][0]

with xr.open_dataset(f'{MESH}/fesom.mesh.diag.nc', decode_times=False) as m:
    mlat = m['lat'].values
    if np.abs(mlat).max() < 4:
        mlat = np.rad2deg(mlat)
    sa = m['nod_area'].values[0, :]

def model_mld(var):
    acc, n = None, 0
    for y in YEARS:
        f = glob.glob(f'{root}/outdata/fesom/{var}.fesom.{y}.nc')
        if not f:
            continue
        with xr.open_dataset(f[0], decode_times=False) as d:
            nm = var if var in d.data_vars else list(d.data_vars)[-1]
            a = np.abs(d[nm].values[FALL, :]).mean(axis=0)
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n) if n else None

with xr.open_dataset(WOA, decode_times=False) as d:
    wmld = np.squeeze(d['M_an'].values)
    wlat, wlon = d['lat'].values, d['lon'].values
print(f'arm {ARM}, model months Oct-Dec of {YEARS[0]}-{YEARS[-1]}')
print("reference: WOA18 ocean_mixed_layer_thickness Fall (Oct-Dec) 1981-2010, 1 deg\n")

m1, m2 = model_mld('MLD1'), model_mld('MLD2')
hdr = (f'  {"band":24s} {"WOA18":>9s} {"MLD1":>9s} {"MLD1-obs":>10s} '
       f'{"MLD2":>9s} {"MLD2-obs":>10s}')
print(hdr); print('  ' + '-' * (len(hdr) - 2))
for name, la0, la1 in BANDS:
    ws = (wlat >= la0) & (wlat < la1)
    sub = wmld[ws, :]
    w = np.broadcast_to(np.cos(np.deg2rad(wlat[ws]))[:, None], sub.shape)
    ok = np.isfinite(sub)
    wo = float(np.average(sub[ok], weights=w[ok])) if ok.any() else np.nan
    ms = (mlat >= la0) & (mlat < la1)
    row = [wo]
    for mm in (m1, m2):
        if mm is None:
            row += [np.nan, np.nan]; continue
        v = float(np.nansum(mm[ms] * sa[ms]) / np.nansum(sa[ms]))
        row += [v, v - wo]
    print(f'  {name:24s} {row[0]:9.1f} {row[1]:9.1f} {row[2]:+10.1f} '
          f'{row[3]:9.1f} {row[4]:+10.1f}')
