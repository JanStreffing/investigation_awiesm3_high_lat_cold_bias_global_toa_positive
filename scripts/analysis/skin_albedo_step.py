"""Why the ice-skin commits cost 1.6 K of Arctic winter: the melting-albedo step.

THE HYPOTHESIS.  ice_thermo_cpl.F90:700 selects the ice/snow albedo with a hard step at
the freezing point -- t < 273.15 gives albsn/albi, t >= 273.15 gives albsnm/albim -- and
FESOM ships the result to OpenIFS as A_Ice_albedo (gen_forcing_couple.F90:389). The skin
is capped at 273.15 immediately after the solve, so "melting" means "sitting exactly on
the cap". If 11X's self-referential solve parked on that cap more often than 13A's
OIFS-anchored skin does, 11X ran with a spuriously LOW ice albedo, and the fix removes
that -- which is the observed 0.470 -> 0.528 jump, not anything to do with winter.

WHAT IS MEASURED.  For each arm and each hemisphere's melt season, over ice-covered nodes
only, the fraction of node-days with the skin on the melting branch. Also the cold tail,
because the same broken solve produced both shoulders.

TRAPS.  ist and a_ice are DAILY (365/366); m_ice and the OIFS fields are MONTHLY. Indexing
a daily array with a month number silently returns the first twelve days of January and
has already reversed one conclusion in this campaign.

Usage:  python3 scripts/analysis/skin_albedo_step.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
ARMS=['11X','13C','13D','13A','11Y']
YEARS=(1357,1358,1359)
AICE_MIN=0.15          # ice-covered: matches the extent definition used elsewhere

def doy(m0,m1,nt):
    ml=np.array([31,28,31,30,31,30,31,31,30,31,30,31]);  ml[1]=29 if nt==366 else 28
    e=np.cumsum(ml); s=np.r_[0,e[:-1]]
    return np.r_[tuple(np.arange(s[m],e[m]) for m in range(m0,m1))]

print(f'melting-branch occupancy (skin >= 273.14 K) over nodes with a_ice >= {AICE_MIN}')
print(f'years {YEARS[0]}-{YEARS[-1]}\n')
hdr=f"{'arm':>5} | {'NH JJA melt%':>12} {'NH JJA mean':>11} | {'SH DJF melt%':>12} {'SH DJF mean':>11} | {'NH cold<220K%':>13}"
print(hdr); print('-'*len(hdr))
for arm in ARMS:
    acc={k:[] for k in ('nm','nt','sm','st','nc')}
    for y in YEARS:
        pi,pa=f'{R}/{arm}/outdata/fesom/ist.fesom.{y}.nc', f'{R}/{arm}/outdata/fesom/a_ice.fesom.{y}.nc'
        if not (os.path.exists(pi) and os.path.exists(pa)): break
        with xr.open_dataset(pi,decode_times=False) as di, xr.open_dataset(pa,decode_times=False) as da:
            t=np.asarray(di['ist'].values,dtype=np.float64)
            a=np.asarray(da['a_ice'].values,dtype=np.float64)
            lat=np.asarray(di['lat'].values)
            nt=t.shape[0]
            for hemi,mons,key in (('N',(5,8),'n'),('S',(11,12),'s')):
                idx = doy(*mons,nt) if hemi=='N' else np.r_[doy(11,12,nt),doy(0,2,nt)]
                m   = (lat>=60) if hemi=='N' else (lat<=-60)
                tt,aa = t[np.ix_(idx,m)], a[np.ix_(idx,m)]
                k = (aa>=AICE_MIN) & np.isfinite(tt)
                if k.any():
                    acc[key+'m'].append(float((tt[k]>=273.14).mean()))
                    acc[key+'t'].append(float(tt[k].mean()))
            idx=np.r_[doy(11,12,nt),doy(0,2,nt)]; m=lat>=60
            tt,aa=t[np.ix_(idx,m)],a[np.ix_(idx,m)]
            k=(aa>=AICE_MIN)&np.isfinite(tt)
            if k.any(): acc['nc'].append(float((tt[k]<220.0).mean()))
    if not acc['nm']: print(f'{arm:>5} | missing'); continue
    f=lambda k: np.mean(acc[k])
    print(f'{arm:>5} | {100*f("nm"):11.2f}% {f("nt"):11.2f} | {100*f("sm"):11.2f}% {f("st"):11.2f} | {100*f("nc"):12.3f}%')

# ---------------------------------------------------------------------------
# Part 2: the melt-pond amplifier.
# albsn->albsnm is only 0.81->0.77, too small on its own for the observed
# 0.470->0.528 box-albedo jump. But use_meltponds=.true., and meltpond_area
# (ice_thermo_cpl.F90:512) is driven by the top-melt rate and by `t` through the
# lid freeze/melt, so a skin that parks on the melting cap grows ponds. Pond
# albedo is ~0.28 against snow's 0.81, so pond fraction is the real lever.
#
# TRAP: apnd/hpnd are MONTHLY (12) while a_ice is DAILY. Aggregate the daily
# field to months rather than indexing either one with the other's time axis.
_ML=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
def to_monthly(a):
    nt=a.shape[0]
    if nt==12: return a
    ml=_ML.copy();  ml[1]=29 if nt==366 else 28
    out=np.empty((12,)+a.shape[1:]); i=0
    for m in range(12):
        out[m]=np.nanmean(a[i:i+ml[m]],axis=0); i+=ml[m]
    return out

print()
print('melt ponds, NH JJA (months 6-8), nodes with a_ice >= 0.15')
h2=f"{'arm':>5} | {'apnd mean':>10} {'apnd|pond':>10} {'nodes ponded':>13} {'hpnd|pond m':>12}"
print(h2); print('-'*len(h2))
JJA=[5,6,7]
for arm in ARMS:
    A_,Ac,Af,H_=[],[],[],[]
    for y in YEARS:
        pp,pa,ph=(f'{R}/{arm}/outdata/fesom/{v}.fesom.{y}.nc' for v in ('apnd','a_ice','hpnd'))
        if not all(os.path.exists(q) for q in (pp,pa,ph)): break
        with xr.open_dataset(pp,decode_times=False) as dp, xr.open_dataset(pa,decode_times=False) as da, \
             xr.open_dataset(ph,decode_times=False) as dh:
            ap=to_monthly(np.asarray(dp['apnd'].values,dtype=np.float64))
            hp=to_monthly(np.asarray(dh['hpnd'].values,dtype=np.float64))
            ai=to_monthly(np.asarray(da['a_ice'].values,dtype=np.float64))
            m=np.asarray(dp['lat'].values)>=60
            ap,hp,ai=ap[np.ix_(JJA,m)],hp[np.ix_(JJA,m)],ai[np.ix_(JJA,m)]
            k=(ai>=AICE_MIN)&np.isfinite(ap)
            if k.any():
                A_.append(float(ap[k].mean())); Af.append(float((ap[k]>0.01).mean()))
                pk=k&(ap>0.01)
                Ac.append(float(ap[pk].mean()) if pk.any() else 0.0)
                H_.append(float(hp[pk].mean()) if pk.any() else 0.0)
    if not A_: print(f'{arm:>5} | missing'); continue
    print(f'{arm:>5} | {np.mean(A_):10.4f} {np.mean(Ac):10.4f} {100*np.mean(Af):12.2f}% {np.mean(H_):12.4f}')
