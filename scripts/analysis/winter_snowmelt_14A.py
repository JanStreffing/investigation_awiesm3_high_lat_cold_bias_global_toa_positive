"""Where and at what skin temperature does 14A melt snow in Feb-Apr?

With snowmelt_tgate=.false. 14A melts NH snow at -0.08 to -0.09 m/month in Feb-Apr 1350.
If that happens in the sunless high Arctic at skins far below 0 C it is an artefact of
testing melt with a flux the atmosphere evaluated at a cold skin; if it sits at the
sunlit marginal ice zone it may be real.  Split by latitude band and report the skin
temperature at the nodes doing the melting.

TRAPS.  ist, a_ice DAILY; thdgrsnw MONTHLY.

Usage:  python3 scripts/analysis/winter_snowmelt_14A.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
_ML=np.array([31,28,31,30,31,30,31,31,30,31,30,31]); SEC=_ML*86400.0
MN=['Jan','Feb','Mar','Apr','May','Jun']
def path(arm,v):
    for p in (f'{R}/{arm}/outdata/fesom/{v}.fesom.1350.nc',
              f'{R}/{arm}/run_13500101-13591231/work/{v}.fesom_1350-1350.nc'):
        if os.path.exists(p): return p
def monthly(a):
    if a.shape[0]==12: return a
    e=np.cumsum(_ML); s=np.r_[0,e[:-1]]; k=int((e<=a.shape[0]).sum())
    return np.stack([np.nanmean(a[s[m]:e[m]],0) for m in range(k)])
def rd(arm,v):
    with xr.open_dataset(path(arm,v),decode_times=False) as d:
        return np.asarray(d[v].values,dtype=np.float64), np.asarray(d['lat'].values)
for arm in ('11Y','15A','15B','15C'):
    ts,lat=rd(arm,'thdgrsnw'); ist,_=rd(arm,'ist'); ai,_=rd(arm,'a_ice')
    ts=monthly(ts); ist=monthly(ist); ai=monthly(ai)
    print(f'\n{arm} 1350: snow melt (m/month, band mean over all nodes) and skin where melting')
    for lo,hi in ((60,70),(70,80),(80,90)):
        m=(lat>=lo)&(lat<hi); w=np.cos(np.deg2rad(lat[m]))
        row=[];sk=[]
        for k in range(6):
            x=ts[k][m]; f=np.isfinite(x)
            row.append(np.average(x[f],weights=w[f])*SEC[k] if f.any() else np.nan)
            mel=f&(x*SEC[k]<-0.01)&(ai[k][m]>=0.15)&np.isfinite(ist[k][m])
            sk.append(np.average(ist[k][m][mel],weights=w[mel]) if mel.any() else np.nan)
        print(f'  {lo}-{hi}N melt  '+' '.join(f'{MN[k]} {row[k]:+.3f}' for k in range(6)))
        print(f'  {lo}-{hi}N skin  '+' '.join(f'{MN[k]} {sk[k]:6.1f}' for k in range(6))
              +'   K, at nodes melting >1 cm/month')
