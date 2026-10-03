"""Southern Ocean state and ocean heat uptake, fast part, 4-year window across arms.

For each arm, 1350-1353 (reads outdata, or a running leg's work dir):
  * global ocean heat uptake from thetaoga: rho*cp*V*dtheta/dt, W/m2 of Earth area
    (annual means, least-squares slope over the window), to set against net TOA
  * surface heat ENTRY by band from fh (stored POSITIVE UPWARD, Round 35), W/m2 of band
  * Southern Ocean winter (JAS) MLD2: band means and convective area (MLD2 > 1000 m)
  * SH sea-ice volume and mean floe thickness in September (m_ice, a_ice; ratio of means)
The depth-resolved heat content and PHC3 comparison are in so_heat_depth_4yr.py.

Usage:  ARMS=11X,11Y,13A,15A,15B,15C python3 scripts/analysis/so_heat_fast_4yr.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESH='/work/ab0246/a270092/input/fesom2/core3/mesh.nc'   # CORE3 geometry (lat, cell_area, depth_bnds)
ARMS=os.environ.get('ARMS','11X,11Y,13A,15A,15B,15C').split(',')
YEARS=range(int(os.environ.get('Y0',1350)),int(os.environ.get('Y1',1353))+1)
RHO,CP,AE=1025.0,3992.0,5.101e14
_ML=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

with xr.open_dataset(MESH,decode_times=False) as m:
    lat=m['lat'].values.astype('f8'); A0=m['cell_area'].values.astype('f8')
    db=np.abs(m['depth_bnds'].values.astype('f8')); nlv=m['depth_lev'].values.astype(int)
VOL=float((A0*db[np.clip(nlv,0,len(db)-1)]).sum())     # cell area x bottom depth
print(f'CORE3: {lat.size} nodes, ocean volume {VOL:.4e} m3, surface {A0.sum():.4e} m2')

def path(arm,v,y):
    for p in (f'{R}/{arm}/outdata/fesom/{v}.fesom.{y}.nc',
              f'{R}/{arm}/run_13500101-13591231/work/{v}.fesom_{y}-{y}.nc'):
        if os.path.exists(p): return p
def rd(arm,v,y):
    p=path(arm,v,y)
    if p is None: return None
    with xr.open_dataset(p,decode_times=False) as d:
        return np.asarray(d[v].values,dtype=np.float64)
def monthly(a):
    n=a.shape[0]
    if n==12: return a
    ml=_ML.copy(); ml[1]=29 if n==366 else 28; e=np.cumsum(ml); s=np.r_[0,e[:-1]]
    return np.stack([np.nanmean(a[s[k]:e[k]],0) for k in range(12)])
bm=lambda x,sel: float(np.nansum(x[sel]*A0[sel])/A0[sel][np.isfinite(x[sel])].sum())

BANDS=[('60-90N',60,90),('30-60N',30,60),('30S-30N',-30,30),('45-30S',-45,-30),('60-45S',-60,-45),('90-60S',-90,-60)]
res={}
for arm in ARMS:
    r={}
    th=[np.nanmean(np.squeeze(x)) for x in (rd(arm,'thetaoga',y) for y in YEARS) if x is not None]
    if len(th)>=2:
        slope=np.polyfit(np.arange(len(th)),th,1)[0]           # K / yr
        r['OHU thetaoga']=RHO*CP*VOL*slope/(365.25*86400)/AE
    fh=[np.nanmean(x,0) for x in (rd(arm,'fh',y) for y in YEARS) if x is not None]
    if fh:
        f=-np.mean(fh,0)                                      # positive INTO the ocean
        r['fh global (Earth)']=float(np.nansum(f*A0))/AE
        for n,lo,hi in BANDS: r[f'fh {n}']=bm(f,(lat>=lo)&(lat<hi))
    mld=[monthly(x)[[6,7,8]].mean(0) for x in (rd(arm,'MLD2',y) for y in YEARS) if x is not None]
    if mld:
        M=np.abs(np.mean(mld,0))
        for n,lo,hi in (('60-45S',-60,-45),('90-60S',-90,-60)): r[f'MLD2 JAS {n}']=bm(M,(lat>=lo)&(lat<hi))
        so=lat<-50
        r['conv area >1000m']=float(A0[so&(M>1000)].sum())/1e12
        r['conv area >2000m']=float(A0[so&(M>2000)].sum())/1e12
    vs=[];ts=[]
    for y in YEARS:
        mi,ai=rd(arm,'m_ice',y),rd(arm,'a_ice',y)
        if mi is None or ai is None: continue
        mi,ai=monthly(mi)[8],monthly(ai)[8]; sh=lat<0
        vs.append(float(np.nansum(mi[sh]*A0[sh]))/1e12)
        ts.append(float(np.nansum(mi[sh]*A0[sh])/np.nansum(ai[sh]*A0[sh])))
    if vs: r['SH ice vol Sep 1e3km3']=np.mean(vs); r['SH floe thk Sep m']=np.mean(ts)
    res[arm]=r
keys=list(dict.fromkeys(k for r in res.values() for k in r))
print(f'\nwindow {YEARS[0]}-{YEARS[-1]}; W/m2 positive INTO ocean; MLD m; area M km2')
print(f"{'metric':<24}"+''.join(f'{a:>10}' for a in ARMS))
for k in keys: print(f'{k:<24}'+''.join(f"{res[a].get(k,np.nan):10.3f}" for a in ARMS))
