# Recovered from the Claude Code session transcript on 2026-07-28.
# Original location was the ephemeral session scratchpad, which was wiped.
# Session c56eada6-56b5-4b15-83cd-c6ed69cf48d9, written 07-28T07:31:51.

import numpy as np, xarray as xr, glob, json, warnings; warnings.filterwarnings('ignore')
SP='/tmp/claude-24456/-work-ab0246-a270092-postprocessing-investigation-awiesm3-high-lat-cold-bias-global-toa-positive/c56eada6-56b5-4b15-83cd-c6ed69cf48d9/'
RT='/work/bb1469/a270270/runtime/awiesm3-v3.4/'; OURS='/work/bb1469/a270092/runtime/awiesm3-v3.4/'
BOXES={'Boreal 55-70N':(55,70,-180,180),'Siberia':(55,75,60,180),'E. Siberia':(55,75,90,160)}
LSM=xr.open_dataset(RT+'Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')['lsm'].isel(time_counter=0)
def boxmean(da,lsm,box):
    la0,la1,lo0,lo1=box
    lat=da.lat.values; lon=((da.lon.values+180)%360)-180
    yi=(lat>=la0)&(lat<=la1); xi=(lon>=lo0)&(lon<=lo1)
    v=da.values[:,yi,:][:,:,xi]; L=lsm[yi,:][:,xi]>0.5
    w=np.cos(np.deg2rad(lat[yi]))[:,None]*np.ones(L.shape); w=np.where(L,w,0.0)
    return (v*w).sum(axis=(1,2))/w.sum()
out=json.load(open(SP+'scratchpad/seasonal.json'))
# new AMIP baseline (correct build): regular grid, yrs 1872-1879 (drop 2 spin-up yrs)
A='/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_base/outdata/oifs/'
a2=xr.open_dataset(A+'atm_remapped_1m_2t_1m_1875-1875.nc')['2t']
assert a2.shape[1:]==LSM.shape, (a2.shape,LSM.shape)
acc={b:np.zeros(12) for b in BOXES}; n=0
for y in range(1872,1880):
    ds=xr.open_dataset(A+f'atm_remapped_1m_2t_1m_{y}-{y}.nc')
    for b,box in BOXES.items(): acc[b]+=boxmean(ds['2t'],LSM.values,box)
    ds.close(); n+=1
out['AMIP baseline (new build)']={b:(acc[b]/n-273.15).tolist() for b in BOXES}
print('AMIP new build done',flush=True)
for lab,path in [('09A = baseline + newSeaIce',RT+'Tuning_test_09A_lpjguess_Baseline_coupled_fromCRUNCEP_newSeaIce'),
                 ('09B = 06T + newSeaIce',RT+'Tuning_test_09B_06T_1hCPL_MOSPP_KPPLOW_CRUNCEPinit_newSeaIce'),
                 ('09C = 06V + newSeaIce',OURS+'Tuning_test_09C_06V_CRUNCEPinit_newSeaIce')]:
    acc={b:np.zeros(12) for b in BOXES}; n=0
    for y in range(1370,1380):
        f=f'{path}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc'
        ds=xr.open_dataset(f)
        for b,box in BOXES.items(): acc[b]+=boxmean(ds['2t'],LSM.values,box)
        ds.close(); n+=1
    out[lab]={b:(acc[b]/n-273.15).tolist() for b in BOXES}; print(lab,'done',flush=True)
json.dump(out,open(SP+'scratchpad/seasonal.json','w'),indent=1)
ref=out['CRUNCEP3 reference']
for s in ['AMIP (atmosphere only)','AMIP baseline (new build)','080a (CRUNCEP veg init)','09A = baseline + newSeaIce','08B = 06V','09C = 06V + newSeaIce','09B = 06T + newSeaIce']:
    if s not in out: continue
    print(f"{s:30s} "+" | ".join(f"{b[:11]} JJA {np.array(out[s][b])[5:8].mean()-np.array(ref[b])[5:8].mean():+.2f} ANN {np.mean(out[s][b])-np.mean(ref[b]):+.2f}" for b in BOXES))
