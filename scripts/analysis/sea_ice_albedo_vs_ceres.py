"""How bright is the model's sea ice, against CERES?

DEFINITION MATCHED ON BOTH SIDES.  All-sky broadband SURFACE albedo over consolidated sea
ice, flux-weighted:  alpha = sum(SW_up * w) / sum(SW_down * w)  over cells whose OWN ice
product says concentration >= CI_MIN, in the given band and months.
  * model:  SW_down = ssrd, SW_up = ssrd - ssr (ssr is NET), mask on OIFS `ci`
  * obs:    sfc_sw_up_all_clim / sfc_sw_down_all_clim (CERES EBAF Ed4.1), mask on HadISST
            sea-ice concentration climatology (same 1x1 grid as CERES)
Flux weighting is essential: averaging the per-cell ratio blows up as SW_down -> 0 near the
polar night terminator and produces a number with no physical meaning.

Masking each side by its own ice concentration is deliberate. It answers "where each has
consolidated ice, how bright is that ice", which is the surface-property question, and it
does NOT confound the answer with the model's 44% extent error.

CAVEATS THAT MUST TRAVEL WITH THE NUMBER.
  * CERES EBAF-Surface is a radiative-transfer product constrained by TOA, not a direct
    surface measurement. Over bright, cloudy, low-sun polar scenes it is among its weakest
    regimes. Treat differences under ~0.03 as within its uncertainty.
  * CERES/HadISST are present-day; the model is 1850 preindustrial. Present-day Arctic ice
    is younger, thinner and more ponded, so the observed albedo is if anything biased LOW
    relative to a true preindustrial. That works AGAINST finding the model too bright.
  * `fal` is reported alongside as the model's own diagnosed grid-box albedo, a cross-check
    on the flux ratio, not an independent measurement.

Usage:  python3 scripts/analysis/sea_ice_albedo_vs_ceres.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R   = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CER = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
SIC = '/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc'
ARMS=['11X','13C','13D','13A','11Y']
YEARS=(1357,1358,1359)
CI_MIN=0.90
SEASONS=[('NH Apr-Sep',[3,4,5,6,7,8],60,90),('NH Jun-Aug',[5,6,7],60,90),
         ('SH Oct-Mar',[9,10,11,0,1,2],-90,-60),('SH Dec-Feb',[11,0,1],-90,-60)]

def fluxalb(up,dn,ci,lat2d,lo,hi):
    w=np.cos(np.deg2rad(lat2d))
    k=(ci>=CI_MIN)&(lat2d>=lo)&(lat2d<hi)&np.isfinite(up)&np.isfinite(dn)&(dn>0)
    if not k.any(): return np.nan,0.0
    return float((up[k]*w[k]).sum()/(dn[k]*w[k]).sum()), float(w[k].sum())

# ---- observations -----------------------------------------------------------
with xr.open_dataset(CER,decode_times=False) as dc, xr.open_dataset(SIC,decode_times=False) as ds:
    clat=np.squeeze(dc['lat'].values); clon=np.squeeze(dc['lon'].values)
    cup=np.squeeze(dc['sfc_sw_up_all_clim'].values); cdn=np.squeeze(dc['sfc_sw_down_all_clim'].values)
    slat=np.squeeze(ds['latitude'].values); sic=np.squeeze(ds['sic'].values)
    if np.nanmax(sic)>1.5: sic=sic/100.0
    if slat[0]>slat[-1]:                       # HadISST runs N->S, CERES S->N
        sic=sic[:,::-1,:]; slat=slat[::-1]
    assert abs(slat[0]-clat[0])<1.0, 'obs grids not aligned in latitude'
    clat2=np.broadcast_to(clat[:,None],cup.shape[1:])
    # HadISST longitudes start at -179.5, CERES at 0.5: roll to match
    sic=np.roll(sic,180,axis=2)
    obs={}
    for name,mons,lo,hi in SEASONS:
        num=den=0.0
        w=np.cos(np.deg2rad(clat2))
        for m in mons:
            k=(sic[m]>=CI_MIN)&(clat2>=lo)&(clat2<hi)&np.isfinite(cup[m])&(cdn[m]>0)
            if k.any(): num+=float((cup[m][k]*w[k]).sum()); den+=float((cdn[m][k]*w[k]).sum())
        obs[name]=num/den if den>0 else np.nan

# ---- model ------------------------------------------------------------------
def mload(arm,var,y):
    p=f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(p): return None,None
    with xr.open_dataset(p,decode_times=False) as d:
        k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values), np.squeeze(d['lat'].values)

print(f'All-sky surface albedo over ice with concentration >= {CI_MIN:.2f}, flux-weighted')
print(f'model years {YEARS[0]}-{YEARS[-1]}; obs = CERES EBAF Ed4.1 surface / HadISST SIC\n')
hdr=f"{'season':>12} {'CERES':>7} "+' '.join(f'{a:>7}' for a in ARMS)+f"  {'13A-obs':>8}"
print(hdr); print('-'*len(hdr))
res={}
for name,mons,lo,hi in SEASONS:
    row=[]
    for arm in ARMS:
        vals=[]
        for y in YEARS:
            ssr,lat=mload(arm,'ssr',y); ssrd,_=mload(arm,'ssrd',y); ci,_=mload(arm,'ci',y)
            if ssr is None or ssrd is None or ci is None: break
            lat2=np.broadcast_to(lat[:,None],ssr.shape[1:])
            num=den=0.0; w=np.cos(np.deg2rad(lat2))
            for m in mons:
                up=ssrd[m]-ssr[m]; dn=ssrd[m]
                k=(ci[m]>=CI_MIN)&(lat2>=lo)&(lat2<hi)&np.isfinite(up)&(dn>0)
                if k.any(): num+=float((up[k]*w[k]).sum()); den+=float((dn[k]*w[k]).sum())
            if den>0: vals.append(num/den)
        row.append(np.mean(vals) if vals else np.nan)
    res[name]=row
    d13=row[ARMS.index('13A')]-obs[name]
    print(f"{name:>12} {obs[name]:7.3f} "+' '.join(f'{v:7.3f}' for v in row)+f"  {d13:+8.3f}")

# ---------------------------------------------------------------------------
# Hardening.  Snow-covered land also sits near 0.8, so a coastal cell leaking into
# the mask would manufacture exactly the bias reported above. Re-run with an explicit
# ocean mask, report the mean concentration inside each mask (open water at 0.06 drags
# a grid-box albedo down, so the masks must be comparable), and cross-check the flux
# ratio against the model's own diagnosed grid-box albedo `fal`.
LSMF=('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
with xr.open_dataset(LSMF,decode_times=False) as d:
    lsm=np.squeeze(d['lsm'].values); lsm=lsm[0] if lsm.ndim==3 else lsm
ocean=lsm<=0.5
print('\nocean-masked, with mask composition and the model\'s own `fal`:')
hdr2=(f"{'season':>12} {'CERES':>7} {'obs ci':>7} | "
      +' '.join(f'{a:>7}' for a in ARMS)+f" | {'13A ci':>7} {'13A fal':>8}")
print(hdr2); print('-'*len(hdr2))
for name,mons,lo,hi in SEASONS:
    # obs mask composition
    w=np.cos(np.deg2rad(clat2)); sn=sd=0.0
    for m in mons:
        k=(sic[m]>=CI_MIN)&(clat2>=lo)&(clat2<hi)&np.isfinite(cup[m])&(cdn[m]>0)
        if k.any(): sn+=float((sic[m][k]*w[k]).sum()); sd+=float(w[k].sum())
    obs_ci=sn/sd if sd>0 else np.nan
    row=[]; ci13=fal13=np.nan
    for arm in ARMS:
        vals=[]; cis=[]; fals=[]
        for y in YEARS:
            ssr,lat=mload(arm,'ssr',y); ssrd,_=mload(arm,'ssrd',y)
            ci,_=mload(arm,'ci',y);     fal,_=mload(arm,'fal',y)
            if ssr is None or ssrd is None or ci is None: break
            lat2=np.broadcast_to(lat[:,None],ssr.shape[1:]); w2=np.cos(np.deg2rad(lat2))
            num=den=0.0; cn=cd=0.0; fn=0.0
            for m in mons:
                up=ssrd[m]-ssr[m]; dn=ssrd[m]
                k=(ci[m]>=CI_MIN)&ocean&(lat2>=lo)&(lat2<hi)&np.isfinite(up)&(dn>0)
                if k.any():
                    num+=float((up[k]*w2[k]).sum()); den+=float((dn[k]*w2[k]).sum())
                    cn+=float((ci[m][k]*w2[k]).sum()); cd+=float(w2[k].sum())
                    if fal is not None: fn+=float((fal[m][k]*w2[k]).sum())
            if den>0:
                vals.append(num/den); cis.append(cn/cd)
                if fal is not None: fals.append(fn/cd)
        row.append(np.mean(vals) if vals else np.nan)
        if arm=='13A':
            ci13=np.mean(cis) if cis else np.nan
            fal13=np.mean(fals) if fals else np.nan
    print(f"{name:>12} {obs[name]:7.3f} {obs_ci:7.3f} | "
          +' '.join(f'{v:7.3f}' for v in row)+f" | {ci13:7.3f} {fal13:8.3f}")

# ---------------------------------------------------------------------------
# Monthly shape, and the energy the gap is worth.
print('\nNH 60-90N monthly, flux-weighted albedo over ci>=0.90 ocean, + SW down and')
print('the absorbed-SW difference the albedo gap implies (model minus CERES, W/m2):')
mn=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
print(f"{'':>10} "+' '.join(f'{m:>6}' for m in mn[2:10]))
def nh_month(arm):
    A=[];D=[]
    for m in range(12):
        va=[];vd=[]
        for y in YEARS:
            ssr,lat=mload(arm,'ssr',y); ssrd,_=mload(arm,'ssrd',y); ci,_=mload(arm,'ci',y)
            if ssr is None: break
            lat2=np.broadcast_to(lat[:,None],ssr.shape[1:]); w2=np.cos(np.deg2rad(lat2))
            up=ssrd[m]-ssr[m]; dn=ssrd[m]
            k=(ci[m]>=CI_MIN)&ocean&(lat2>=60)&np.isfinite(up)&(dn>0)
            if k.any():
                va.append(float((up[k]*w2[k]).sum()/(dn[k]*w2[k]).sum()))
                vd.append(float((dn[k]*w2[k]).sum()/w2[k].sum())/3600.0)
        A.append(np.mean(va) if va else np.nan); D.append(np.mean(vd) if vd else np.nan)
    return np.array(A),np.array(D)
oA=[];oD=[]
w=np.cos(np.deg2rad(clat2))
for m in range(12):
    k=(sic[m]>=CI_MIN)&(clat2>=60)&np.isfinite(cup[m])&(cdn[m]>0)
    oA.append(float((cup[m][k]*w[k]).sum()/(cdn[m][k]*w[k]).sum()) if k.any() else np.nan)
    oD.append(float((cdn[m][k]*w[k]).sum()/w[k].sum()) if k.any() else np.nan)
oA=np.array(oA); oD=np.array(oD)
print(f"{'CERES':>10} "+' '.join(f'{v:6.3f}' for v in oA[2:10]))
print(f"{'SWdn obs':>10} "+' '.join(f'{v:6.1f}' for v in oD[2:10]))
for arm in ('11X','13A'):
    A,D=nh_month(arm)
    print(f"{arm:>10} "+' '.join(f'{v:6.3f}' for v in A[2:10]))
    print(f"{'  SWdn':>10} "+' '.join(f'{v:6.1f}' for v in D[2:10]))
    dq=-(A-oA)*D
    print(f"{'  dSWabs':>10} "+' '.join(f'{v:6.1f}' for v in dq[2:10])
          +f"   | Jun-Aug mean {np.nanmean(dq[5:8]):6.1f} W/m2")

# ---------------------------------------------------------------------------
# Is it too bright BECAUSE there is too much ice, or too bright at the same ice?
# The ci>=0.90 masks above match concentration but not LOCATION: the model's
# consolidated pack reaches further south into stronger sun and into cells the
# observations say are open water. Co-locate to separate the two.
#   BOTH      model ci>=0.90 AND obs sic>=0.90   -> intrinsic surface brightness
#   MODEL ONLY  model consolidated, obs not      -> the extent error's own albedo
#   OBS ONLY    obs consolidated, model not
# Model is nearest-neighbour binned from its 192x400 grid onto the CERES 1x1.
print('\nco-located decomposition, NH 60-90N, Jun-Aug, flux-weighted')
_, mlat = mload('13A','ci',YEARS[0])
with xr.open_dataset(f'{R}/13A/outdata/oifs/atm_remapped_1m_ci_{YEARS[0]}-{YEARS[0]}.nc',
                     decode_times=False) as d:
    mlon=np.squeeze(d['lon'].values)
iy=np.abs(clat[:,None]-mlat[None,:]).argmin(axis=1)
ix=np.abs(((clon[:,None]-mlon[None,:]+180)%360-180)).argmin(axis=1)
rg=lambda a: a[np.ix_(iy,ix)]
oc1=rg(ocean)
hdr3=f"{'arm':>5} {'class':>11} {'area %':>7} {'model a':>8} {'CERES a':>8} {'diff':>7} {'SWdn':>6} {'dSWabs':>7}"
print(hdr3); print('-'*len(hdr3))
for arm in ('11X','13A'):
    tot=None; store={}
    for cls in ('BOTH','MODEL ONLY','OBS ONLY'):
        mn_=md_=on_=od_=aw=0.0
        for y in YEARS:
            ssr,_=mload(arm,'ssr',y); ssrd,_=mload(arm,'ssrd',y); ci,_=mload(arm,'ci',y)
            if ssr is None: break
            for m in (5,6,7):
                mci=rg(ci[m]); up=rg(ssrd[m]-ssr[m]); dn=rg(ssrd[m])
                oci=sic[m]; oup=cup[m]; odn=cdn[m]
                w=np.cos(np.deg2rad(clat2))
                base=(clat2>=60)&oc1&np.isfinite(up)&(dn>0)&np.isfinite(oup)&(odn>0)
                if   cls=='BOTH':       k=base&(mci>=CI_MIN)&(oci>=CI_MIN)
                elif cls=='MODEL ONLY': k=base&(mci>=CI_MIN)&(oci< CI_MIN)
                else:                   k=base&(mci< CI_MIN)&(oci>=CI_MIN)
                if k.any():
                    mn_+=float((up[k]*w[k]).sum()); md_+=float((dn[k]*w[k]).sum())
                    on_+=float((oup[k]*w[k]).sum()); od_+=float((odn[k]*w[k]).sum())
                    aw+=float(w[k].sum())
        store[cls]=(mn_/md_ if md_>0 else np.nan, on_/od_ if od_>0 else np.nan,
                    md_/aw/3600.0 if aw>0 else np.nan, aw)
    tot=sum(v[3] for v in store.values())
    for cls,(ma,oa,swd,aw) in store.items():
        print(f"{arm:>5} {cls:>11} {100*aw/tot:6.1f}% {ma:8.3f} {oa:8.3f} {ma-oa:+7.3f} "
              f"{swd:6.1f} {-(ma-oa)*swd:+7.1f}")
