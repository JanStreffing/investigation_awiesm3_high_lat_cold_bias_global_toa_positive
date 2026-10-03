"""First-summer check of the 14-series: does the snow on Arctic sea ice melt now?

The 14 arms remove the coupled snow-melt gate (and add an albedo ramp / ECHAM albedos).
The mechanism acts within the first melt season, so this reads the running work
directories for 1350 and compares against 11Y's 1350 (same commit, same precision).
Only months fully written are used: the files are still being appended.

Quantities, NH 60-90N, ice nodes with a_ice >= 0.15:
  snow depth on the floe  m_snow/a_ice   (FESOM albedo leaves the snow branch below 1 mm)
  melt pond fraction      apnd
TRAPS.  a_ice is DAILY, m_snow and apnd MONTHLY; aggregate before combining.

Usage:  python3 scripts/analysis/early_snow_check_14.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
Y=1350
_ML=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
MN=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

def path(arm,v):
    # running legs write XIOS names ({v}.fesom_Y-Y.nc); renamed at leg end
    for p in (f'{R}/{arm}/run_13500101-13591231/work/{v}.fesom_{Y}-{Y}.nc',
              f'{R}/{arm}/run_13500101-13591231/work/{v}.fesom.{Y}.nc',
              f'{R}/{arm}/outdata/fesom/{v}.fesom.{Y}.nc'):
        if os.path.exists(p): return p
    return None

def monthly_daily(a):
    """Daily -> monthly over the COMPLETE months present."""
    n=a.shape[0]; e=np.cumsum(_ML); s=np.r_[0,e[:-1]]
    k=int((e<=n).sum())
    return np.stack([np.nanmean(a[s[m]:e[m]],0) for m in range(k)]) if k else None

for arm in ('11Y','14A','14B','14C'):
    pa,ps,pp=path(arm,'a_ice'),path(arm,'m_snow'),path(arm,'apnd')
    if not all((pa,ps,pp)): print(f'{arm}: files missing'); continue
    with xr.open_dataset(pa,decode_times=False) as da, xr.open_dataset(ps,decode_times=False) as ds, \
         xr.open_dataset(pp,decode_times=False) as dp:
        ai=monthly_daily(np.asarray(da['a_ice'].values,dtype=np.float64))
        ms=np.asarray(ds['m_snow'].values,dtype=np.float64)
        ap=np.asarray(dp['apnd'].values,dtype=np.float64)
        lat=np.asarray(da['lat'].values)
    if ai is None: print(f'{arm}: no complete month yet'); continue
    nm=min(ai.shape[0],ms.shape[0],ap.shape[0])
    m=lat>=60
    sn=[];pn=[]
    for k in range(nm):
        sel=(ai[k][m]>=0.15)&np.isfinite(ms[k][m])
        sn.append(float((ms[k][m][sel]/np.maximum(ai[k][m][sel],1e-6)).mean()) if sel.any() else np.nan)
        pn.append(float(np.nanmean(ap[k][m][sel])) if sel.any() else np.nan)   # apnd is NaN off-ice
    print(f'\n{arm}  ({nm} complete months of {Y})')
    print('   month   '+' '.join(f'{MN[k]:>6}' for k in range(3,nm)))
    print('   snow m  '+' '.join(f'{sn[k]:6.3f}' for k in range(3,nm)))
    print('   ponds   '+' '.join(f'{pn[k]:6.3f}' for k in range(3,nm)))

# ---------------------------------------------------------------------------
# Part 2: the snow BUDGET, without the ratio.  m_snow/a_ice over a_ice>=0.15 nodes rises
# spuriously when floe area melts out from under a snow pack that is not removed with it,
# and the node selection changes month to month.  Here: band-mean grid quantities over ALL
# nodes 60-90N (cos-lat weighted), plus FESOM's thermodynamic snow change (thdgrsnw, <0 =
# melt) and the snowfall reaching the ice/ocean (snow), both as m/month.
print('\n=== snow budget, NH 60-90N band means over all nodes (no ice selection) ===')
SEC=_ML*86400.0
def band(arm,v,daily=False):
    p=path(arm,v)
    if p is None: return None
    with xr.open_dataset(p,decode_times=False) as d:
        if v not in d.data_vars: return None
        a=np.asarray(d[v].values,dtype=np.float64); lat=np.asarray(d['lat'].values)
    if a.shape[0]>12: a=monthly_daily(a)
    m=lat>=60; w=np.cos(np.deg2rad(lat[m]))
    out=[]
    for k in range(a.shape[0]):
        x=a[k][m]; f=np.isfinite(x)
        out.append(float(np.average(x[f],weights=w[f])) if f.any() else np.nan)
    return np.array(out)
for arm in ('11Y','14A','14B','14C'):
    ai,ms=band(arm,'a_ice'),band(arm,'m_snow')
    ts,sf=band(arm,'thdgrsnw'),band(arm,'snow')
    if ai is None or ms is None: print(f'{arm}: missing'); continue
    n=min(len(ai),len(ms)); r=range(3,min(n,10))
    print(f'\n{arm}   '+' '.join(f'{MN[k]:>7}' for k in r))
    print('  a_ice     '+' '.join(f'{ai[k]:7.3f}' for k in r))
    print('  m_snow m  '+' '.join(f'{ms[k]:7.4f}' for k in r))
    if ts is not None:
        print('  thdgrsnw  '+' '.join(f'{ts[k]*SEC[k]:+7.4f}' for k in r if k<len(ts))+'   m/month, <0 melt')
    if sf is not None:
        print('  snowfall  '+' '.join(f'{sf[k]*SEC[k]:+7.4f}' for k in r if k<len(sf))+'   m/month (FESOM snow flux)')
