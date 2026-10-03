"""Co-located JJA sea-ice albedo vs CERES for the first summer (1350) of the 14 arms.

Same definition as sea_ice_albedo_vs_ceres.py (flux-weighted all-sky surface albedo,
cells where BOTH model ci >= 0.90 and HadISST sic >= 0.90, ocean only, model
nearest-neighbour binned onto the CERES 1x1).  Reads running work dirs for 14A/B/C and
outdata for 11Y (same commit and precision, the control).  One year is noisy for climate
but the albedo is a direct surface property and responds within the season.

Usage:  python3 scripts/analysis/early_albedo_check_14.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
CER='/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
SIC='/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc'
LSMF=('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
CI=0.90; Y=1350; JJA=(5,6,7)

with xr.open_dataset(CER,decode_times=False) as dc, xr.open_dataset(SIC,decode_times=False) as ds:
    clat=np.squeeze(dc['lat'].values); clon=np.squeeze(dc['lon'].values)
    cup=np.squeeze(dc['sfc_sw_up_all_clim'].values); cdn=np.squeeze(dc['sfc_sw_down_all_clim'].values)
    slat=np.squeeze(ds['latitude'].values); sic=np.squeeze(ds['sic'].values)
if np.nanmax(sic)>1.5: sic=sic/100.0
if slat[0]>slat[-1]: sic=sic[:,::-1,:]
sic=np.roll(sic,180,axis=2)
clat2=np.broadcast_to(clat[:,None],cup.shape[1:]); W=np.cos(np.deg2rad(clat2))
with xr.open_dataset(LSMF,decode_times=False) as d:
    lsm=np.squeeze(d['lsm'].values); lsm=lsm[0] if lsm.ndim==3 else lsm
    mlat=np.squeeze(d['lat'].values); mlon=np.squeeze(d['lon'].values)
iy=np.abs(clat[:,None]-mlat[None,:]).argmin(1)
ix=np.abs(((clon[:,None]-mlon[None,:]+180)%360-180)).argmin(1)
rg=lambda a: a[np.ix_(iy,ix)]; oc=rg(lsm<=0.5)

def mload(arm,v):
    for p in (f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{v}_{Y}-{Y}.nc',
              f'{R}/{arm}/run_13500101-13591231/work/atm_remapped_1m_{v}_{Y}-{Y}.nc'):
        if os.path.exists(p):
            with xr.open_dataset(p,decode_times=False) as d:
                k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
                return np.squeeze(d[k].values)
    return None

print(f'co-located JJA {Y}, NH 60-90N, model ci>={CI} AND obs sic>={CI}, ocean, flux-weighted')
print(f"{'arm':>5} {'model a':>8} {'CERES a':>8} {'diff':>7} {'dSWabs W/m2':>12} {'cells':>7}")
for arm in ('11Y','14A','15A','15B','15C'):
    ssr,ssrd,ci=mload(arm,'ssr'),mload(arm,'ssrd'),mload(arm,'ci')
    if ssr is None or ssrd is None or ci is None or ssr.shape[0]<8:
        print(f'{arm:>5}   not yet available'); continue
    mn=md=on=od=aw=0.0
    for m in JJA:
        up=rg(ssrd[m]-ssr[m]); dn=rg(ssrd[m]); mc=rg(ci[m])
        k=(clat2>=60)&oc&np.isfinite(up)&(dn>0)&(mc>=CI)&(sic[m]>=CI)&(cdn[m]>0)
        mn+=(up[k]*W[k]).sum(); md+=(dn[k]*W[k]).sum()
        on+=(cup[m][k]*W[k]).sum(); od+=(cdn[m][k]*W[k]).sum(); aw+=W[k].sum()
    ma,oa=mn/md,on/od; sw=md/aw/3600.0
    print(f'{arm:>5} {ma:8.3f} {oa:8.3f} {ma-oa:+7.3f} {-(ma-oa)*sw:+12.1f} {aw:7.0f}')
