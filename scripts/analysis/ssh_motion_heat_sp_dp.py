"""Is the SP surface heat term real free-surface motion, or a W inconsistent with it?

Under linfs the tracer's top vertical velocity should equal the free-surface tendency,
so the heat the surface term carries can be rebuilt from the model's own daily ssh:
    H = -vcpw * sum( deta/dt * T_top * area ) / A_earth      (W/m2 of Earth area)
with deta/dt from the day-1 to day-5 daily means and T_top the matching mean sst.
The heat probe measured the surface term (Iadvv) at SP -0.51 and DP -0.01 (5-day means,
HP_SP2 / HP_DP2, same restart).  If H reproduces that, the SP model really moves its free
surface differently; if H is the same in both, the SP tracer W is inconsistent with the
free surface it should follow.  Split by region (cavity nodes from mesh.nc cav_nod_mask).

Usage:  python3 scripts/analysis/ssh_motion_heat_sp_dp.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings, glob
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'; VCPW=4.2e6; AE=5.101e14
with xr.open_dataset('/work/ab0246/a270092/input/fesom2/core3/mesh.nc',decode_times=False) as m:
    A=m['cell_area'].values.astype('f8'); lat=m['lat'].values.astype('f8')
    cav=(m['cav_nod_mask'].values>0) if 'cav_nod_mask' in m else np.zeros_like(A,bool)
reg={'cavity':cav,'open 90-60S':(~cav)&(lat<-60),'open 60S-60N':(~cav)&(lat>=-60)&(lat<60),'open 60-90N':(~cav)&(lat>=60)}
print(f'cavity nodes: {cav.sum()} of {cav.size}')
def rd(e,v):
    f=glob.glob(f'{R}/{e}/run_*/work/{v}.fesom_1350-1350.nc')[0]
    with xr.open_dataset(f,decode_times=False) as d:
        return np.asarray(d[v].values,dtype='f8')
res={}
for e in ('HP_SP2','HP_DP2'):
    eta=rd(e,'ssh'); sst=rd(e,'sst')
    detadt=(eta[4]-eta[0])/(4*86400.0); T=0.5*(sst[0]+sst[4])
    k=np.isfinite(detadt)&np.isfinite(T)
    res[e]=(detadt,T)
    tot=-VCPW*np.nansum((detadt*T*A)[k])/AE
    vol=np.nansum((detadt*A)[k])
    parts=' '.join(f"{n}={-VCPW*np.nansum((detadt*T*A)[k&msk])/AE:+.3f}" for n,msk in reg.items())
    print(f'{e}: H={tot:+.3f} W/m2  (sum deta/dt*area = {vol:+.3e} m3/s)   {parts}')
(dS,TS),(dD,TD)=res['HP_SP2'],res['HP_DP2']
k=np.isfinite(dS)&np.isfinite(dD)
dd=dS-dD
print(f'SP-DP deta/dt: rms {np.sqrt(np.nanmean(dd[k]**2)):.3e} m/s, area-mean {np.nansum((dd*A)[k])/A[k].sum():+.3e} m/s')
print('SP-DP heat by region: '+' '.join(f"{n}={-VCPW*np.nansum((dd*TS*A)[k&msk])/AE:+.3f}" for n,msk in reg.items()))
