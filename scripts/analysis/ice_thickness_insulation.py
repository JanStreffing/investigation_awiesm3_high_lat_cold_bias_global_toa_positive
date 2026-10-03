"""Does the too-thick ice throttle the winter ocean-to-atmosphere heat flux?

THE HYPOTHESIS (operator, 2026-09-10).  Conduction through sea ice goes as
k*(T_freeze - T_skin)/h, so ice that is too thick insulates the ocean.  In winter that
means (a) too little heat reaching the atmosphere, which IS the DJF high-latitude cold
bias, and (b) heat retained in the ocean, which is the accumulation the heat-routing work
measured going below 300 m.  If it holds, the campaign's two objectives are one mechanism.

WHAT IS MEASURED, per arm, 60-90N and 90-60S, over nodes with a_ice >= 0.15:
  * ice thickness m_ice/a_ice (true floe thickness, not the grid-cell mean)
  * qcon, the conductive flux through the ice
  * fh, the total surface heat flux (POSITIVE UPWARD in this configuration -- verified in
    the Round 35 free checks; getting this sign wrong inverts the conclusion)
  * sst under the ice, as the reservoir the conduction is drawing from

TRAPS.  a_ice and ist are DAILY (365/366); m_ice, qcon, fh, sst are MONTHLY.  m_ice is
volume per unit grid area, so dividing by a_ice is required to get floe thickness --
comparing m_ice against an observed thickness compares different quantities.

Usage:  python3 scripts/analysis/ice_thickness_insulation.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
ARMS=['11X','13C','13D','13A','11Y']
YEARS=(1357,1358,1359)
AMIN=0.15
_ML=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

def to_monthly(a):
    nt=a.shape[0]
    if nt==12: return a
    ml=_ML.copy(); ml[1]=29 if nt==366 else 28
    o=np.empty((12,)+a.shape[1:]); i=0
    for m in range(12): o[m]=np.nanmean(a[i:i+ml[m]],axis=0); i+=ml[m]
    return o

def get(arm,y,var):
    p=f'{R}/{arm}/outdata/fesom/{var}.fesom.{y}.nc'
    if not os.path.exists(p): return None,None
    with xr.open_dataset(p,decode_times=False) as d:
        k=var if var in d.data_vars else [c for c in d.data_vars if 'bnds' not in c][0]
        return to_monthly(np.asarray(d[k].values,dtype=np.float64)), np.asarray(d['lat'].values)

for hemi,latsel,mons,lbl in (('NH',lambda l:l>=60,[0,1,11],'DJF'),
                             ('SH',lambda l:l<=-60,[5,6,7],'JJA')):
    print(f'\n=== {hemi} {lbl}, nodes with a_ice >= {AMIN} ===')
    hdr=(f"{'arm':>5} {'thick m':>8} {'a_ice':>7} {'m_ice m':>8} {'qcon W/m2':>10} "
         f"{'fh W/m2':>9} {'sst C':>7} {'k dT/h est':>11}")
    print(hdr); print('-'*len(hdr))
    for arm in ARMS:
        acc={k:[] for k in ('th','ai','mi','qc','fh','sst','est')}
        for y in YEARS:
            ai,lat=get(arm,y,'a_ice')
            if ai is None: break
            mi,_=get(arm,y,'m_ice'); qc,_=get(arm,y,'qcon')
            fh,_=get(arm,y,'fh');    ss,_=get(arm,y,'sst')
            ist,_=get(arm,y,'ist')
            m=latsel(lat); w=np.cos(np.deg2rad(lat[m]))
            for k_ in mons:
                a=ai[k_,m]; sel=(a>=AMIN)&np.isfinite(a)
                if not sel.any(): continue
                ww=w[sel]
                th=np.where(a>1e-6, mi[k_,m]/np.maximum(a,1e-6), np.nan)
                acc['th'].append(np.average(th[sel][np.isfinite(th[sel])],
                                            weights=ww[np.isfinite(th[sel])]))
                acc['ai'].append(np.average(a[sel],weights=ww))
                acc['mi'].append(np.average(mi[k_,m][sel],weights=ww))
                for key,arr in (('qc',qc),('fh',fh),('sst',ss)):
                    if arr is None: continue
                    v_=arr[k_,m][sel]; f_=np.isfinite(v_)
                    if f_.any(): acc[key].append(np.average(v_[f_],weights=ww[f_]))
                if ist is not None:
                    t_=ist[k_,m][sel]; h_=th[sel]
                    g=np.isfinite(t_)&np.isfinite(h_)&(h_>0.05)
                    if g.any():
                        acc['est'].append(np.average(2.1656*(271.35-t_[g])/h_[g],
                                                     weights=ww[g]))
        if not acc['th']: print(f'{arm:>5}  missing'); continue
        f=lambda k: np.mean(acc[k]) if acc[k] else np.nan
        print(f"{arm:>5} {f('th'):8.2f} {f('ai'):7.3f} {f('mi'):8.2f} {f('qc'):10.2f} "
              f"{f('fh'):9.2f} {f('sst'):7.3f} {f('est'):11.1f}")
