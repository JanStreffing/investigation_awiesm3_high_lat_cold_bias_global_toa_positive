"""Snow on NH sea ice as a BAND budget, not a per-floe ratio.

Why this exists: mean(m_snow/a_ice) over nodes with a_ice >= 0.15 is dominated by
low-concentration nodes, where a snow pack left on vanishing ice divides by a small
concentration. It reported 13A "carrying 0.27-0.32 m of snow through July-August"; band
means show 11Y's snow melting in June. Use ratio-of-means (band m_snow / band a_ice) and
the grid quantities themselves.

Usage:  ARMS=11X:1357-1359,13A:1357-1359 python3 scripts/analysis/snow_budget_band.py
        (default: 11X and 13A 1357-59, plus 11Y and 14A 1350 for Jan-Mar)
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
_ML=np.array([31,28,31,30,31,30,31,31,30,31,30,31]); SEC=_ML*86400.0
MN=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
SPEC=os.environ.get('ARMS','11X:1357-1359,13A:1357-1359,11Y:1350-1350,14A:1350-1350')

def path(arm,v,y):
    for p in (f'{R}/{arm}/outdata/fesom/{v}.fesom.{y}.nc',
              f'{R}/{arm}/run_13500101-13591231/work/{v}.fesom_{y}-{y}.nc'):
        if os.path.exists(p): return p
    return None
def monthly(a):
    n=a.shape[0]
    if n==12: return a
    ml=_ML.copy(); ml[1]=29 if n==366 else 28
    e=np.cumsum(ml); s=np.r_[0,e[:-1]]; k=int((e<=n).sum())
    return np.stack([np.nanmean(a[s[m]:e[m]],0) for m in range(k)])
def band(arm,v,years):
    rows=[]
    for y in years:
        p=path(arm,v,y)
        if p is None: continue
        with xr.open_dataset(p,decode_times=False) as d:
            if v not in d.data_vars: return None
            a=monthly(np.asarray(d[v].values,dtype=np.float64)); lat=np.asarray(d['lat'].values)
        m=lat>=60; w=np.cos(np.deg2rad(lat[m])); r=np.full(12,np.nan)
        for k in range(a.shape[0]):
            x=a[k][m]; f=np.isfinite(x)
            if f.any(): r[k]=np.average(x[f],weights=w[f])
        rows.append(r)
    return np.nanmean(np.array(rows),0) if rows else None

for item in SPEC.split(','):
    arm,yr=item.split(':'); y0,y1=map(int,yr.split('-')); ys=range(y0,y1+1)
    ai,ms,ts=band(arm,'a_ice',ys),band(arm,'m_snow',ys),band(arm,'thdgrsnw',ys)
    if ai is None or ms is None: print(f'\n{arm} {yr}: missing'); continue
    print(f'\n{arm} {yr}, NH 60-90N band means (all nodes)')
    print('               '+' '.join(f'{m:>7}' for m in MN[:9]))
    print('  a_ice        '+' '.join(f'{x:7.3f}' for x in ai[:9]))
    print('  m_snow m     '+' '.join(f'{x:7.4f}' for x in ms[:9]))
    print('  snow/floe m  '+' '.join(f'{x:7.4f}' for x in (ms/ai)[:9])+'   ratio of means')
    if ts is not None:
        print('  thdgrsnw     '+' '.join(f'{x*s:+7.4f}' for x,s in zip(ts[:9],SEC[:9]))+'   m/month, <0 melt')
