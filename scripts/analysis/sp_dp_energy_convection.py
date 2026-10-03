"""Does single precision leak ocean heat, and does it switch off Southern Ocean convection?

The 4-year scorecard (so_heat_fast_4yr.py) split cleanly by FESOM precision: the SP arms
(11Y, 15A-C) lose 0.6-0.7 W/m2 more ocean heat (thetaoga) than their surface flux (fh)
delivers, the DP arms (11X, 13A) close to within 0.1-0.2; and only the DP arms convect
in the Southern Ocean.  Four years cannot separate a systematic effect from timing and
noise, so this goes year by year over the long runs:
  closure  = ocean heat uptake from thetaoga (per-year difference of annual means)
             minus the global surface flux into the ocean (-fh, annual mean)   [W/m2 Earth]
  conv     = winter (JAS mean) area south of 50S with MLD2 > 1000 m            [M km2]
11X is DP at 4d7da170, 11Y SP at 0b07f9a3, 13A DP at 0b07f9a3 (13A vs 11Y = precision only).

Usage:  python3 scripts/analysis/sp_dp_energy_convection.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
RHO,CP,AE,SY=1025.0,3992.0,5.101e14,365.25*86400
with xr.open_dataset('/work/ab0246/a270092/input/fesom2/core3/mesh.nc',decode_times=False) as m:
    lat=m['lat'].values.astype('f8'); A0=m['cell_area'].values.astype('f8')
    db=np.abs(m['depth_bnds'].values.astype('f8')); nlv=m['depth_lev'].values.astype(int)
VOL=float((A0*db[np.clip(nlv,0,len(db)-1)]).sum()); so=lat<-50
def rd(arm,v,y):
    # outdata, or the work dir of a running or cancelled leg (XIOS naming)
    for p in (f'{R}/{arm}/outdata/fesom/{v}.fesom.{y}.nc',
              f'{R}/{arm}/run_13500101-13591231/work/{v}.fesom_{y}-{y}.nc'):
        if os.path.exists(p):
            with xr.open_dataset(p,decode_times=False) as d: a=np.asarray(d[v].values,dtype='f8')
            return a if a.shape[0]>0 else None   # XIOS file of an unfinished year has no records
    return None
ARMS=os.environ.get('ARMS','11X:1350:1390,11Y:1350:1390,13A:1350:1360,15C:1350:1360,15F:1350:1360')
for arm,yrs in [(a,range(int(b),int(c))) for a,b,c in (x.split(':') for x in ARMS.split(','))]:
    th=[];fh=[];cv=[]
    for y in yrs:
        t,f,m=rd(arm,'thetaoga',y),rd(arm,'fh',y),rd(arm,'MLD2',y)
        if t is None or f is None or m is None: break
        th.append(np.nanmean(t)); fh.append(float(np.nansum(-np.nanmean(f,0)*A0))/AE)
        M=np.abs(m[[6,7,8]].mean(0)); cv.append(float(A0[so&(M>1000)].sum())/1e12)
    th=np.array(th); fh=np.array(fh); cv=np.array(cv); n=len(th)
    ohu=RHO*CP*VOL*np.diff(th)/SY/AE                      # between consecutive annual means
    fm=0.5*(fh[1:]+fh[:-1])                               # flux over the same interval
    clo=ohu-fm
    print(f'\n{arm} ({yrs[0]}-{yrs[0]+n-1}, {n} yr)')
    print(f'  mean OHU {ohu.mean():+.3f}  mean fh-in {fm.mean():+.3f}  closure {clo.mean():+.3f} '
          f'(sd {clo.std(ddof=1):.3f}, se {clo.std(ddof=1)/np.sqrt(len(clo)):.3f}) W/m2')
    for a,b in ((0,4),(0,10),(10,20),(20,30),(30,40)):
        if b-1<=len(clo): print(f'    years {a+1:>2}-{b:<2}: closure {clo[a:b-1].mean():+.3f}   conv>1000m {cv[a:b].mean():5.2f} M km2')
    print('  conv area per year: '+' '.join(f'{x:.1f}' for x in cv))
