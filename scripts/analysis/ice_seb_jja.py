"""JJA surface energy balance over consolidated Arctic sea ice, model vs CERES.

QUESTION.  With the snow-melt gate on, snow melts only when OpenIFS's ice skin reaches
melting, and in 13A that happens in ~15 % of JJA node-days while the observed pack sits
at 0 C all summer.  Is the surface losing energy in summer, and through which term?

DEFINITION.  Cells where the model has ci >= 0.90 AND HadISST sic >= 0.90 (ocean), model
nearest-neighbour binned onto the CERES 1x1, cos-lat weighted, JJA.  IFS fluxes are
ACCUMULATED J/m2 per output hour (/3600) and positive DOWNWARD.  CERES gives SW and LW
only; turbulent terms are model-only.  CERES is present-day against a 1850 model, so the
LW-down comparison carries a forcing offset of a few W/m2.

Usage:  RUNS=11X:1357-1359,13A:1357-1359,11Y:1350-1350,14A:1350-1350 python3 scripts/analysis/ice_seb_jja.py
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
CI=0.90; JJA=(5,6,7)
RUNS=os.environ.get('RUNS','11X:1357-1359,13A:1357-1359,11Y:1350-1350,14A:1350-1350')

with xr.open_dataset(CER,decode_times=False) as dc, xr.open_dataset(SIC,decode_times=False) as ds:
    clat=np.squeeze(dc['lat'].values); clon=np.squeeze(dc['lon'].values)
    C={k:np.squeeze(dc[v].values) for k,v in (('swd','sfc_sw_down_all_clim'),('swu','sfc_sw_up_all_clim'),
                                              ('lwd','sfc_lw_down_all_clim'),('lwu','sfc_lw_up_all_clim'))}
    slat=np.squeeze(ds['latitude'].values); sic=np.squeeze(ds['sic'].values)
if np.nanmax(sic)>1.5: sic=sic/100.0
if slat[0]>slat[-1]: sic=sic[:,::-1,:]
sic=np.roll(sic,180,axis=2)
clat2=np.broadcast_to(clat[:,None],sic.shape[1:]); W=np.cos(np.deg2rad(clat2))
with xr.open_dataset(LSMF,decode_times=False) as d:
    lsm=np.squeeze(d['lsm'].values); lsm=lsm[0] if lsm.ndim==3 else lsm
    mlat=np.squeeze(d['lat'].values); mlon=np.squeeze(d['lon'].values)
iy=np.abs(clat[:,None]-mlat[None,:]).argmin(1)
ix=np.abs(((clon[:,None]-mlon[None,:]+180)%360-180)).argmin(1)
rg=lambda a: a[np.ix_(iy,ix)]; oc=rg(lsm<=0.5)

def mload(arm,v,y):
    for p in (f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc',
              f'{R}/{arm}/run_13500101-13591231/work/atm_remapped_1m_{v}_{y}-{y}.nc'):
        if os.path.exists(p):
            with xr.open_dataset(p,decode_times=False) as d:
                k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
                return np.squeeze(d[k].values)
    return None

TERMS=('SWdn','SWnet','LWdn','LWnet','SH','LH','NET')
print('JJA, consolidated NH ice (model ci>=0.9 AND obs sic>=0.9), W/m2, positive DOWN')
print(f"{'run':>14} "+' '.join(f'{t:>7}' for t in TERMS)+f" {'albedo':>7} {'cells':>6}")
for item in RUNS.split(','):
    arm,yr=item.split(':'); y0,y1=map(int,yr.split('-'))
    acc={t:0.0 for t in TERMS}; obs={t:0.0 for t in ('SWdn','SWnet','LWdn','LWnet')}; aw=0.0; ok=True
    for y in range(y0,y1+1):
        F={v:mload(arm,v,y) for v in ('ssrd','ssr','strd','str','sshf','slhf','ci')}
        if any(F[v] is None for v in F) or F['ssr'].shape[0]<8: ok=False; break
        for m in JJA:
            mc=rg(F['ci'][m])
            k=(clat2>=60)&oc&(mc>=CI)&(sic[m]>=CI)&np.isfinite(rg(F['ssr'][m]))&(C['swd'][m]>0)
            w=W[k]; aw+=w.sum()
            f=lambda v: (rg(F[v][m])[k]/3600.0*w).sum()
            acc['SWdn']+=f('ssrd'); acc['SWnet']+=f('ssr'); acc['LWdn']+=f('strd'); acc['LWnet']+=f('str')
            acc['SH']+=f('sshf'); acc['LH']+=f('slhf')
            obs['SWdn']+=(C['swd'][m][k]*w).sum(); obs['SWnet']+=((C['swd'][m]-C['swu'][m])[k]*w).sum()
            obs['LWdn']+=(C['lwd'][m][k]*w).sum(); obs['LWnet']+=((C['lwd'][m]-C['lwu'][m])[k]*w).sum()
    if not ok or aw==0: print(f'{arm+" "+yr:>14}  not available'); continue
    for t in acc: acc[t]/=aw
    for t in obs: obs[t]/=aw
    acc['NET']=acc['SWnet']+acc['LWnet']+acc['SH']+acc['LH']
    alb=1-acc['SWnet']/acc['SWdn']
    print(f"{arm+' '+yr:>14} "+' '.join(f'{acc[t]:7.1f}' for t in TERMS)+f" {alb:7.3f} {aw:6.0f}")
    print(f"{'  CERES same':>14} "+' '.join(f'{obs[t]:7.1f}' if t in obs else f"{'':>7}" for t in TERMS)
          +f" {1-obs['SWnet']/obs['SWdn']:7.3f}")
