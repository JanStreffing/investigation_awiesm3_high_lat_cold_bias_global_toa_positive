"""Where the heat goes in the first years, by depth and band, and the Southern Ocean vs PHC3.

Heat-content change between the 1350 and the last-year annual-mean temperature fields,
rho*cp*sum(dT * nod_area[k] * dz[k]), per latitude band and depth range, as W/m2 of
EARTH area (so rows add to the global ocean uptake and compare with net TOA).  All arms
start from the same 1350-01-01 restart, so the difference is what each arm did.
Then the last-year Southern Ocean zonal-mean temperature at fixed depths against PHC3,
binned on each mesh's own latitudes (no interpolation; same 47-level vertical grid) --
method from so_subsurface_vs_phc3.py.  PHC3 is a late-20th-century climatology: the
comparison carries an epoch offset of a few tenths of a degree in the upper ocean.

Usage:  ARMS=11X,11Y,13A,15A,15B,15C Y1=1353 python3 scripts/analysis/so_heat_depth_4yr.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESH='/work/ab0246/a270092/input/fesom2/core3/mesh.nc'   # CORE3 geometry; no per-level area, so
# below-bottom points are masked where temp is NaN or exactly 0 (FESOM fill)
PHC='/work/ab0246/a270092/obs/phc3'
ARMS=os.environ.get('ARMS','11X,11Y,13A,15A,15B,15C').split(',')
Y0,Y1=1350,int(os.environ.get('Y1',1353))
RHO,CP,AE,SY=1025.0,3992.0,5.101e14,365.25*86400
BANDS=[('60-90N',60,90),('30-60N',30,60),('30S-30N',-30,30),('45-30S',-45,-30),('60-45S',-60,-45),('90-60S',-90,-60)]
DEPTHS=[('0-100',0,100),('100-700',100,700),('700-2000',700,2000),('>2000',2000,1e9)]
SOB=[('60-45S',-60,-45),('90-60S',-90,-60)]; SOD=[0,100,500,1000,2000]

def meta(diag):
    with xr.open_dataset(diag,decode_times=False) as d:
        lat=np.asarray(d['lat'].values) if 'lat' in d else np.asarray(d['nodes'].values)[1]
        if np.nanmax(np.abs(lat))<3.2: lat=np.rad2deg(lat)
        return lat,np.asarray(d['nod_area'].values),(np.abs(np.asarray(d['nz'].values)) if 'nz' in d else None)
with xr.open_dataset(MESH,decode_times=False) as m:
    lat=m['lat'].values.astype('f8'); carea=m['cell_area'].values.astype('f8')
    db=np.abs(m['depth_bnds'].values.astype('f8'))

def path(arm,y):
    for p in (f'{R}/{arm}/outdata/fesom/temp.fesom.{y}.nc',
              f'{R}/{arm}/run_13500101-13591231/work/temp.fesom_{y}-{y}.nc'):
        if os.path.exists(p): return p
def ann(arm,y):
    p=path(arm,y)
    if p is None: return None
    with xr.open_dataset(p,decode_times=False) as d:
        v=d['temp']; zd=[c for c in v.dims if c.startswith('nz')][0]
        zv=np.abs(d[zd].values.astype('f8')) if zd in d.variables else None
        a=v.mean('time').transpose(zd,...).values.astype(np.float64)   # (nz, nod2)
    a[a==0.0]=np.nan                                                   # FESOM below-bottom fill
    return a,zv

with xr.open_dataset(f'{PHC}/temp.fesom.1958.nc',decode_times=False) as d:
    phc=np.squeeze(d['temp'].values); nzv=np.asarray(d['nz1'].values)
olat,oarea,_=meta(f'{PHC}/fesom.mesh.diag.nc')
def zonal(T_nz_nod,la,ar,zvals):   # T (nz, nod), ar (nz,nod) or (nod,) -> {(band,depth): mean}
    out={}
    for b,lo,hi in SOB:
        s=(la>=lo)&(la<hi)
        for dt in SOD:
            iz=int(np.argmin(np.abs(zvals-dt))); v=T_nz_nod[iz][s]
            w=(ar[min(iz,ar.shape[0]-1)] if ar.ndim==2 else ar)[s]
            k=np.isfinite(v)&(w>0)&(v>-5)&(v<40)
            out[(b,dt)]=float(np.average(v[k],weights=w[k])) if k.any() else np.nan
    return out
OBS=zonal(phc.T,olat,oarea,np.abs(nzv))

print(f'heat-content change {Y0} -> {Y1} (annual means), W/m2 of EARTH area\n')
SO={}
for arm in ARMS:
    a0,a1=ann(arm,Y0),ann(arm,Y1)
    if a0 is None or a1 is None: print(f'{arm}: missing'); continue
    (T0,zv),(T1,_)=a0,a1; nl=T0.shape[0]
    dz=np.diff(db)[:nl]; zm=0.5*(db[:-1]+db[1:])[:nl]; zv=zm if zv is None else zv
    vol=carea[None,:]*dz[:,None]
    dH=RHO*CP*np.where(np.isfinite(T0)&np.isfinite(T1),T1-T0,0.0)*vol/((Y1-Y0)*SY*AE)
    print(f'{arm:<6} {"band":<8}'+''.join(f'{d[0]:>10}' for d in DEPTHS)+f'{"total":>10}')
    col=np.zeros(len(DEPTHS))
    for b,lo,hi in BANDS:
        s=(lat>=lo)&(lat<hi); row=[dH[np.ix_((zm>=d0)&(zm<d1),s)].sum() for _,d0,d1 in DEPTHS]; col+=row
        print(f'{"":<6} {b:<8}'+''.join(f'{v:+10.3f}' for v in row)+f'{sum(row):+10.3f}')
    print(f'{"":<6} {"ALL":<8}'+''.join(f'{v:+10.3f}' for v in col)+f'{col.sum():+10.3f}\n')
    SO[arm]=zonal(T1,lat,carea,zv)
print(f'\nSouthern Ocean zonal-mean temperature {Y1} vs PHC3 [degC] (model - PHC3 in brackets)')
for b,_,_ in SOB:
    print(f'\n--- {b} ---\n{"depth":>6} {"PHC3":>7} '+''.join(f'{a:>15}' for a in SO))
    for dt in SOD:
        print(f'{dt:>6} {OBS[(b,dt)]:7.2f} '+''.join(f'{SO[a][(b,dt)]:7.2f} ({SO[a][(b,dt)]-OBS[(b,dt)]:+5.2f})' for a in SO))
