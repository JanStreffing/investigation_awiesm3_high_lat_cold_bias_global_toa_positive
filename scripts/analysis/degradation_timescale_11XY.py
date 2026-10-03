"""How fast does 11Y diverge from 11X, and therefore how long must an isolation run be?

Both arms branch from the SAME 1350 restart, so year-by-year divergence is directly
attributable to the three things that differ (FESOM library version, dp->sp precision, and
the 2026-09-08 ice-skin fixes). This does not separate them -- it sizes the experiment that
would.

The number wanted is the year at which the difference clears the interannual scatter. If it
is year 1, one-year isolation arms suffice and cost ~20 minutes each; if it is decade 2,
they do not and the isolation is expensive.

IFS fluxes are accumulated J/m2 over the output step; divided by 3600 with the incoming
solar asserted near 340 W/m2 as the guard.
"""
import os, numpy as np, xarray as xr, warnings
for _v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_v,'1')
warnings.filterwarnings('ignore')

R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
LSMF=('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
ACC=3600.0; YEARS=range(1350,1390)
with xr.open_dataset(LSMF, decode_times=False) as d:
    m=np.squeeze(d['lsm'].values); m=m[0] if m.ndim==3 else m
    lat,lon=np.squeeze(d['lat'].values),np.squeeze(d['lon'].values)
land,ocean=m>0.5,m<=0.5
W=np.broadcast_to(np.cos(np.deg2rad(lat))[:,None],m.shape).copy()
bnd=lambda lo,hi: np.broadcast_to(((lat>=lo)&(lat<hi))[:,None],m.shape)
cell=np.broadcast_to((6.371e6**2*np.cos(np.deg2rad(lat))*np.deg2rad(abs(lat[1]-lat[0]))
                      *2*np.pi/m.shape[1])[:,None],m.shape)
def am(f,s):
    k=s&np.isfinite(f); return float(np.average(f[k],weights=W[k])) if k.any() else np.nan
def yr(arm,var,y):
    p=f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p,decode_times=False) as d:
        k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values)

S={}
for arm in ('11X','11Y'):
    rec={k:[] for k in ('yr','nhext','shext','t2m_glob','t2m_arc','toa','alb')}
    for y in YEARS:
        ci=yr(arm,'ci',y); t2=yr(arm,'2t',y)
        tsr=yr(arm,'tsr',y); ttr=yr(arm,'ttr',y); fal=yr(arm,'fal',y)
        if ci is None or t2 is None: continue
        t2=t2-273.15 if np.nanmean(t2)>100 else t2
        rec['yr'].append(y)
        rec['nhext'].append(float(np.nansum(((ci[2]>0.15)*cell)[bnd(0,90)&ocean]))/1e12)
        rec['shext'].append(float(np.nansum(((ci[8]>0.15)*cell)[bnd(-90,0)&ocean]))/1e12)
        rec['t2m_glob'].append(am(t2.mean(0),np.ones_like(m,bool)))
        rec['t2m_arc'].append(am(t2[[11,0,1]].mean(0),bnd(60,90)))
        rec['toa'].append(am((tsr+ttr).mean(0)/ACC,np.ones_like(m,bool)) if tsr is not None else np.nan)
        rec['alb'].append(am(fal.mean(0),bnd(60,90)) if fal is not None else np.nan)
    S[arm]={k:np.array(v) for k,v in rec.items()}

print(__doc__.split('\n')[0]); print()
print(f'{"year":<6}' + ''.join(f'{n:>26}' for n in
      ('NH ice ext Mar','global T2m','DJF T2m 60-90N','sfc albedo 60-90N')))
print(f'{"":<6}' + ''.join(f'{"11X":>8}{"11Y":>9}{"diff":>9}' for _ in range(4)))
print('-'*110)
for i,y in enumerate(S['11X']['yr']):
    if y > 1369 and y % 2: continue
    ln=f'{y:<6}'
    for k in ('nhext','t2m_glob','t2m_arc','alb'):
        a,b=S['11X'][k][i],S['11Y'][k][i]
        ln+=f'{a:>8.2f}{b:>9.2f}{b-a:>+9.2f}'
    print(ln)
print()
for k,nm,unit in (('nhext','NH March ice extent','1e6 km2'),('t2m_glob','global T2m','K'),
                  ('t2m_arc','DJF T2m 60-90N','K'),('alb','sfc albedo 60-90N','-'),
                  ('toa','net TOA','W/m2')):
    a,b=S['11X'][k],S['11Y'][k]
    d=b-a
    sd=np.nanstd(a[:10],ddof=1)                 # 11X interannual scatter, first decade
    thr=1.96*sd*np.sqrt(2.0)                    # single-year paired detection threshold
    first=next((S['11X']['yr'][i] for i in range(len(d))
                if np.isfinite(d[i]) and abs(d[i])>thr), None)
    print(f'{nm:<22} yr1 diff {d[0]:>+8.3f}   thr(1yr) {thr:>6.3f} {unit:<8}'
          f'  first clears: {first}   yr40 diff {d[-1]:>+8.3f}')
