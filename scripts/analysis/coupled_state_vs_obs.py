"""State of the coupled model against observations: temperature, radiation, albedo, ice.

WHAT THIS IS FOR.  The chain table (coupled_chain_1850.py) shows arm-to-arm movement. This
one asks the other question: where does the model actually STAND against observations, so
the campaign objectives have a number attached rather than a direction.

TRAPS OBSERVED (campaign protocol).
  * IFS TOA and surface fluxes are ACCUMULATED J/m2 over the output step. Divided by 3600
    here, and the incoming solar is asserted to land near 340 W/m2 as a guard -- getting
    this wrong makes every flux ~3600x too large and has done so before.
  * Sea ice AREA is not EXTENT. Both are computed and separately labelled; mixing them
    inverted the sign of a conclusion once.
  * ERA5 is an IFS sibling. Usable for T2m and dynamics, NOT an independent arbiter of
    cloud, snow or turbulent flux for this model.
  * CERES is independent, but present-day against 1850 arms: the offset is forcing, not
    model error. Only the arm-to-arm differences are clean.
  * Band means are area-weighted on each grid in its own geometry.

Usage:  ARMS=11X,11Y python3 scripts/analysis/coupled_state_vs_obs.py
"""
import os, glob, sys
for _v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R   = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
E5  = '/work/ab0246/a270092/obs/era5/netcdf'
CER = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
LSMF= ('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
ACC = 3600.0
Y0, Y1 = int(os.environ.get("Y0",1380)), int(os.environ.get("Y1",1389))
DJF, JJA = [11,0,1], [5,6,7]
ARMS = [a.strip() for a in os.environ.get('ARMS','11X,11Y').split(',')]

BANDS=[('60-90N',60,90),('45-60N',45,60),('30-45N',30,45),('30S-30N',-30,30),
       ('45-30S',-45,-30),('60-45S',-60,-45),('90-60S',-90,-60)]

with xr.open_dataset(LSMF, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values); m = m[0] if m.ndim==3 else m
    lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
land, ocean = m>0.5, m<=0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:,None], m.shape).copy()
bnd = lambda lo,hi: np.broadcast_to(((lat>=lo)&(lat<hi))[:,None], m.shape)

def am(f, s):
    k = s & np.isfinite(f)
    return float(np.average(f[k], weights=W[k])) if k.any() else np.nan

def load(arm, var):
    out=[]
    for y in range(Y0,Y1+1):
        p=f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(p):   # a running leg: complete years sit in its work dir
            p=f'{R}/{arm}/run_13500101-13591231/work/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(p): return None
        with xr.open_dataset(p, decode_times=False) as d:
            k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            out.append(np.squeeze(d[k].values))
    return np.stack(out)                       # (years, 12, lat, lon)

# ---------------------------------------------------------------- observations
with xr.open_dataset(f'{E5}/T2M.nc') as d:
    c=d['t2m'].groupby('time.month').mean('time')
    la=[x for x in c.dims if 'lat' in x][0]; lo=[x for x in c.dims if 'lon' in x][0]
    O=np.asarray(c.interp({la:('lat',lat), lo:('lon',np.linspace(0.45,359.55,m.shape[1]))}).values,float)
O = O-273.15 if np.nanmean(O)>100 else O

ceres=None
if os.path.exists(CER):
    with xr.open_dataset(CER, decode_times=False) as d:
        ceres={}
        # the monthly-climatology file names its fields *_clim, the time series *_mon
        for k in ('toa_net_all','toa_sw_all','toa_lw_all','toa_cre_sw','toa_cre_lw'):
            kk=[x for x in (k+'_mon',k+'_clim') if x in d.data_vars]
            if kk:
                a=d[kk[0]]
                cla=[x for x in a.dims if 'lat' in x][0]; clo=[x for x in a.dims if 'lon' in x][0]
                ceres[k]=np.asarray(a.interp({cla:('lat',lat),
                    clo:('lon',np.linspace(0.45,359.55,m.shape[1]))}).values,float)

print(__doc__.split('\n')[0]); print(f'window {Y0}-{Y1}, arms {ARMS}\n')
rows={}
for arm in ARMS:
    t2=load(arm,'2t')
    if t2 is None: print(f'  {arm}: incomplete output, skipped'); continue
    t2=t2-273.15 if np.nanmean(t2)>100 else t2
    tsr=load(arm,'tsr'); ttr=load(arm,'ttr'); tisr=load(arm,'tisr')
    tsrc=load(arm,'tsrc'); ttrc=load(arm,'ttrc'); fal=load(arm,'fal'); ci=load(arm,'ci')
    r={}
    r['T2m DJF 60-90N'] = am(t2[:,DJF].mean((0,1))-O[DJF].mean(0), land&bnd(60,90))
    r['T2m DJF 30-60N'] = am(t2[:,DJF].mean((0,1))-O[DJF].mean(0), land&bnd(30,60))
    r['T2m JJA Siberia']= am(t2[:,JJA].mean((0,1))-O[JJA].mean(0),
                             land&bnd(50,70)&np.broadcast_to(((lon>=60)&(lon<=140))[None,:],m.shape))
    r['T2m ANN land']   = am(t2.mean((0,1))-O.mean(0), land)
    r['T2m ANN ocean']  = am(t2.mean((0,1))-O.mean(0), ocean)
    r['T2m ANN global'] = am(t2.mean((0,1))-O.mean(0), np.ones_like(m,bool))
    r['T2m ANN 30S-30N']= am(t2.mean((0,1))-O.mean(0), bnd(-30,30))
    r['T2m ANN 90-60S'] = am(t2.mean((0,1))-O.mean(0), bnd(-90,-60))
    if tsr is not None:
        net=(tsr+ttr).mean((0,1))/ACC
        sw =(tsr).mean((0,1))/ACC ; lw=(ttr).mean((0,1))/ACC
        inc=(tisr).mean((0,1))/ACC
        assert 320<am(inc,np.ones_like(m,bool))<360, 'incoming solar off -- ACC wrong'
        r['net TOA global']  = am(net, np.ones_like(m,bool))
        r['TOA SW abs glob'] = am(sw,  np.ones_like(m,bool))
        r['TOA LW out glob'] = -am(lw, np.ones_like(m,bool))
        r['SW CRE 45-65S']   = am((tsr-tsrc).mean((0,1))/ACC, bnd(-65,-45))
        r['SW CRE tropics']  = am((tsr-tsrc).mean((0,1))/ACC, bnd(-30,30))
        r['LW CRE tropics']  = am((ttr-ttrc).mean((0,1))/ACC, bnd(-30,30))
        r['planetary albedo']= 1.0-am(sw,np.ones_like(m,bool))/am(inc,np.ones_like(m,bool))
    if fal is not None:
        r['sfc albedo 60-90N'] = am(fal.mean((0,1)), bnd(60,90))
        r['sfc albedo 90-60S'] = am(fal.mean((0,1)), bnd(-90,-60))
    if ci is not None:
        cell=np.broadcast_to((6.371e6**2*np.cos(np.deg2rad(lat))*np.deg2rad(abs(lat[1]-lat[0]))
                              *2*np.pi/m.shape[1])[:,None], m.shape)
        c_=ci.mean(0)
        nh=bnd(0,90)&ocean; sh=bnd(-90,0)&ocean
        r['NH ice AREA Mar']   = float(np.nansum(c_[2][nh]*cell[nh]))/1e12
        r['NH ice EXTENT Mar'] = float(np.nansum(((c_[2]>0.15)*cell)[nh]))/1e12
        r['SH ice AREA Sep']   = float(np.nansum(c_[8][sh]*cell[sh]))/1e12
        r['SH ice EXTENT Sep'] = float(np.nansum(((c_[8]>0.15)*cell)[sh]))/1e12
        r['NH ice EXTENT Sep'] = float(np.nansum(((c_[8]>0.15)*cell)[nh]))/1e12
        r['SH ice EXTENT Feb'] = float(np.nansum(((c_[1]>0.15)*cell)[sh]))/1e12
    rows[arm]=r

keys=list(rows[ARMS[0]].keys()) if rows else []
print(f'{"metric":<22}' + ''.join(f'{a:>13}' for a in rows) +
      (f'{"OBS":>13}' if ceres else ''))
print('-'*(22+13*(len(rows)+(1 if ceres else 0))))
OBSV={}
if ceres:
    G=np.ones_like(m,bool); S0=340.0          # CERES EBAF global-mean incoming solar
    cm=lambda k,msk: am(ceres[k].mean(0),msk) if k in ceres else np.nan
    OBSV['net TOA global']  = cm('toa_net_all',G)
    OBSV['TOA SW abs glob'] = S0-cm('toa_sw_all',G)
    OBSV['TOA LW out glob'] = cm('toa_lw_all',G)
    OBSV['SW CRE 45-65S']   = cm('toa_cre_sw',bnd(-65,-45))
    OBSV['SW CRE tropics']  = cm('toa_cre_sw',bnd(-30,30))
    OBSV['LW CRE tropics']  = cm('toa_cre_lw',bnd(-30,30))
    OBSV['planetary albedo']= cm('toa_sw_all',G)/S0
for k in ('T2m DJF 60-90N','T2m DJF 30-60N','T2m JJA Siberia','T2m ANN land','T2m ANN ocean',
          'T2m ANN global','T2m ANN 30S-30N','T2m ANN 90-60S'): OBSV[k]=0.0   # rows are biases vs ERA5
OBSV['NH ice EXTENT Mar']=15.0; OBSV['SH ice EXTENT Sep']=18.5   # OSI-SAF climatological
OBSV['NH ice EXTENT Sep']=7.0;  OBSV['SH ice EXTENT Feb']=3.0    # satellite era; PI NH Sep ~+2 higher
for k in keys:
    line=f'{k:<22}'+''.join(f'{rows[a][k]:>13.3f}' for a in rows)
    if k in OBSV and np.isfinite(OBSV.get(k,np.nan)): line+=f'{OBSV[k]:>13.2f}'
    print(line)
