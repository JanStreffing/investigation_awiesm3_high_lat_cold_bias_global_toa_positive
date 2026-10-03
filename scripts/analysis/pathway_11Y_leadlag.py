"""How did 11Y's swing happen?  Order of emergence, month by month, 11Y minus 11X.

THE CHAIN UNDER TEST.  If the ice-skin change is the driver, the order must be

    ist rises -> upward longwave over ice rises -> ice grows -> albedo rises
              -> absorbed shortwave falls -> T2m falls

and each link must emerge AFTER the one before it.  If albedo or shortwave moves first,
the driver is not the skin and the chain is wrong.

WHY MONTHLY, AND WHY THE FIRST TWO YEARS.  By year 6 every term has cleared threshold and
they grow together, so the ordering is unrecoverable.  The signal is in the first months,
before the ice-albedo feedback closes and while the terms are still separable.

MIXED OUTPUT FREQUENCIES.  ist and a_ice are written DAILY (365/366 steps), m_ice and the
OIFS fields MONTHLY (12).  Indexing a daily array with a month number silently returns the
first twelve DAYS of January -- it does not fail, it just answers a different question.  All
daily fields are aggregated to monthly means here before anything is combined.

METHOD.  Both arms branch from the SAME 1350 restart on the SAME mesh (core3, 220509), so
a month-by-month difference is attributable without regridding.  Sea-surface fields are
area-weighted over ocean only.  ist is weighted by ice concentration -- an unweighted mean
includes open water where the ice skin temperature means nothing.  IFS fluxes are
accumulated J/m2 over the output step and are divided by 3600.

EMERGENCE is the first month a term exceeds 2 sigma of the control's INTERANNUAL scatter
FOR THAT CALENDAR MONTH.  Taking the scatter of a raw monthly series instead measures the
seasonal cycle, which puts the floor for absorbed shortwave above 100 W/m2 and reports that
nothing ever emerges.  That is a weaker claim than a formal lead-lag correlation and a more
honest one: these series diverge monotonically, so a lagged correlation is dominated by the
shared trend and would manufacture a lead that is not there.

MONTH 1 IS THE CLEAN FINGERPRINT.  Both arms start from an identical restart, so the first
month's difference is the direct, unamplified response to the code change, before chaotic
divergence or any feedback contributes.
"""
import os, numpy as np, xarray as xr, warnings
for _v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_v,'1')
warnings.filterwarnings('ignore')

R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESHF='/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
LSMF=('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
ACC=3600.0; YRS=range(1350,1354)          # 4 years monthly
NOISE=range(1350,1370)                    # 11X years used for the per-calendar-month scatter
m=xr.open_dataset(MESHF,decode_times=False)
flat=m['lat'].values.astype('f8'); farea=m['cell_area'].values.astype('f8')
assert flat.size==220509,'MESH GUARD'
with xr.open_dataset(LSMF,decode_times=False) as d:
    lsm=np.squeeze(d['lsm'].values); lsm=lsm[0] if lsm.ndim==3 else lsm
    alat=np.squeeze(d['lat'].values)
ocean=lsm<=0.5
AW=np.broadcast_to(np.cos(np.deg2rad(alat))[:,None],lsm.shape).copy()
abnd=lambda lo,hi: np.broadcast_to(((alat>=lo)&(alat<hi))[:,None],lsm.shape)

_ML=np.array([31,28,31,30,31,30,31,31,30,31,30,31])

def to_monthly(a):
    """Daily (365/366) -> monthly means. Monthly (12) passes through."""
    nt=a.shape[0]
    if nt==12: return a
    if nt not in (365,366):
        raise ValueError(f'unexpected time length {nt}')
    ml=_ML.copy()
    if nt==366: ml[1]=29
    out=np.empty((12,)+a.shape[1:]); i=0
    for m in range(12):
        out[m]=np.nanmean(a[i:i+ml[m]],axis=0); i+=ml[m]
    return out

def fes(arm,var,y):
    p=f'{R}/{arm}/outdata/fesom/{var}.fesom.{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p,decode_times=False) as d:
        k=[x for x in d.data_vars if x not in('bounds_lon','bounds_lat','time_bounds')][-1]
        return to_monthly(d[k].values.astype('f8'))
def oif(arm,var,y):
    p=f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p,decode_times=False) as d:
        k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values)

def series(arm, lo, hi):
    """monthly series of the chain terms over one band"""
    out={k:[] for k in ('ist','a_ice','m_ice','fal','LWup','SWabs','T2m')}
    fsel=(flat>=lo)&(flat<hi)
    asel=abnd(lo,hi)&ocean
    for y in YRS:
        ist=fes(arm,'ist',y); ai=fes(arm,'a_ice',y); mi=fes(arm,'m_ice',y)
        fal=oif(arm,'fal',y); ssr=oif(arm,'ssr',y); strr=oif(arm,'str',y); t2=oif(arm,'2t',y)
        if ist is None or fal is None: continue
        for mo in range(12):
            w=farea*fsel*np.clip(ai[mo],0,1)              # ice-weighted
            out['ist'].append(np.nansum(ist[mo]*w)/max(np.nansum(w),1e-30))
            wa=farea*fsel
            out['a_ice'].append(np.nansum(ai[mo]*wa)/np.nansum(wa))
            out['m_ice'].append(np.nansum(mi[mo]*wa)/np.nansum(wa))
            k=asel&np.isfinite(fal[mo])
            out['fal'].append(np.average(fal[mo][k],weights=AW[k]))
            out['LWup'].append(-np.average((strr[mo]/ACC)[k],weights=AW[k]))   # +ve = loss
            out['SWabs'].append(np.average((ssr[mo]/ACC)[k],weights=AW[k]))
            out['T2m'].append(np.average(t2[mo][k],weights=AW[k]))
    return {k:np.array(v) for k,v in out.items()}

def series_noise(arm, lo, hi):
    global YRS
    keep=YRS; YRS=NOISE
    try: return series(arm, lo, hi)
    finally: YRS=keep

MON=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
for nm,lo,hi in (('ARCTIC 60-90N',60,90),('ANTARCTIC 90-60S',-90,-60)):
    X=series('11X',lo,hi); Y=series('11Y',lo,hi)
    n=min(len(X['ist']),len(Y['ist']))
    print('='*100); print(f'{nm}   11Y minus 11X, monthly'); print('='*100)
    order=['ist','LWup','m_ice','a_ice','fal','SWabs','T2m']
    lbl ={'ist':'ice skin T [K]','LWup':'LW loss [W/m2]','m_ice':'ice vol [m]',
          'a_ice':'ice frac','fal':'albedo','SWabs':'SW abs [W/m2]','T2m':'T2m [K]'}
    print(f'{"month":<8}'+''.join(f'{lbl[k]:>16}' for k in order))
    d={k:Y[k][:n]-X[k][:n] for k in order}
    for i in range(min(18,n)):
        print(f'{1350+i//12}-{MON[i%12]:<4}'+''.join(f'{d[k][i]:>+16.4f}' for k in order))
    # per-calendar-month interannual scatter of the control
    NX=series_noise('11X',lo,hi)
    print(f'\n{"term":<16}{"month 1":>11}{"emerges":>10}{"thr(2sig)":>11}{"then":>11}{"month 18":>11}')
    for k in order:
        sd=np.array([np.std(NX[k][mo::12],ddof=1) for mo in range(12)])
        first=next((i for i in range(n) if abs(d[k][i])>2*sd[i%12]), None)
        fm=f'{1350+first//12}-{MON[first%12]}' if first is not None else 'never'
        v=f'{d[k][first]:+.4f}' if first is not None else '--'
        print(f'{lbl[k]:<16}{d[k][0]:>+11.4f}{fm:>10}{2*sd[0]:>11.4f}{v:>11}'
              f'{d[k][min(17,n-1)]:>+11.4f}')
    print()
