# Recovered from the Claude Code session transcript on 2026-07-28.
# Original location was the ephemeral session scratchpad, which was wiped.
# Session c56eada6-56b5-4b15-83cd-c6ed69cf48d9, written 07-27T14:33:07.

import numpy as np, xarray as xr, glob, json, warnings
warnings.filterwarnings('ignore')
RT='/work/bb1469/a270270/runtime/awiesm3-v3.4/'
BOXES={'Boreal 55-70N':(55,70,-180,180),'Siberia':(55,75,60,180),'E. Siberia':(55,75,90,160)}
def boxmean_cell(vals, lat, lon, land, box):   # vals (t, cell)
    la0,la1,lo0,lo1=box
    m=(lat>=la0)&(lat<=la1)&(lon>=lo0)&(lon<=lo1)&land
    w=np.cos(np.deg2rad(lat[m]))
    return np.average(vals[:,m],axis=1,weights=w)
def boxmean_grid(da, lsm, box):                # da (t, lat, lon)
    la0,la1,lo0,lo1=box
    lat=da.lat.values; lon=((da.lon.values+180)%360)-180
    yi=(lat>=la0)&(lat<=la1); xi=(lon>=lo0)&(lon<=lo1)
    v=da.values[:,yi,:][:,:,xi]; L=lsm[yi,:][:,xi]>0.5
    w=np.cos(np.deg2rad(lat[yi]))[:,None]*np.ones(L.shape)
    w=np.where(L,w,0.0)
    return (v*w).sum(axis=(1,2))/w.sum()
out={}

# --- cell-grid land mask (shared by AMIP + CRUNCEP3) ---
import cfgrib
g=xr.open_dataset('/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_forcing/input/oifs/95_4/lsm',engine='cfgrib')
land_cell=g['lsm'].values>0.5

# --- AMIP (1870-1879 daily, cell grid) ---
fs=sorted(glob.glob('/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_forcing/outdata/oifs/atm_1d_18*.nc'))
acc={b:np.zeros((12,)) for b in BOXES}; cnt=np.zeros(12)
for f in fs:
    ds=xr.open_dataset(f); lat=ds.lat.values; lon=((ds.lon.values+180)%360)-180
    mon=ds['tas']['time_counter.month'].values; v=ds['tas'].values
    for b,box in BOXES.items():
        s=boxmean_cell(v,lat,lon,land_cell,box)
        for m in range(1,13): acc[b][m-1]+=s[mon==m].mean()
    cnt+=1; ds.close()
out['AMIP (atmosphere only)']={b:(acc[b]/len(fs)-273.15).tolist() for b in BOXES}
print('AMIP done',flush=True)

# --- CRUNCEP3 reference (1901-1910 daily, cell grid) ---
ds=xr.open_dataset('/work/ab0995/a270270/input/cruncep_v7/CRUNCEP_noLPJG_1d_1901-1910_TCO95_calibrated_v3.nc')
lat=ds.lat.values; lon=((ds.lon.values+180)%360)-180
tname=[c for c in ds['tas'].dims if 'time' in c][0]
mon=ds['tas'][tname+'.month'].values if hasattr(ds['tas'],tname) else ds[tname].dt.month.values
ref={}
for b,box in BOXES.items():
    la0,la1,lo0,lo1=box
    m=(lat>=la0)&(lat<=la1)&(lon>=lo0)&(lon<=lo1)&land_cell
    w=np.cos(np.deg2rad(lat[m]))
    sub=ds['tas'][:,m].values
    s=np.average(sub,axis=1,weights=w)
    ref[b]=[float(s[mon==k].mean()-273.15) for k in range(1,13)]
    print('CRUNCEP3',b,flush=True)
ds.close()
out['CRUNCEP3 reference']=ref

# --- coupled runs (1370-1379 monthly, regular grid) ---
RUNS={'06 Baseline':'Tuning_test_06_Baseline',
      '080a (CRUNCEP veg init)':'Tuning_test_080a_lpjguess_Baseline_coupled_fromCRUNCEP',
      '08C = 06T + LPJ lever':'Tuning_test_08C_06T_plus_IBSgreffmin012_CRUNCEPinit',
      '08B = 06V':'Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit',
      '08F = 06V + LPJ levers':'Tuning_test_08F_06V_plus_IBSgreffmin012_plus_C3Gpstemplow12_CRUNCEPinit'}
for lab,d in RUNS.items():
    lsm=xr.open_dataset(f'{RT}{d}/outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')['lsm'].isel(time_counter=0).values
    fs=[f'{RT}{d}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc' for y in range(1370,1380)]
    acc={b:np.zeros(12) for b in BOXES}
    for f in fs:
        ds=xr.open_dataset(f); da=ds['2t']
        for b,box in BOXES.items(): acc[b]+=boxmean_grid(da,lsm,box)
        ds.close()
    out[lab]={b:(acc[b]/len(fs)-273.15).tolist() for b in BOXES}
    print(lab,'done',flush=True)
json.dump(out,open('/tmp/claude-24456/-work-ab0246-a270092-postprocessing-investigation-awiesm3-high-lat-cold-bias-global-toa-positive/c56eada6-56b5-4b15-83cd-c6ed69cf48d9/scratchpad/seasonal.json','w'),indent=1)
print('WROTE')
