# Recovered from the Claude Code session transcript on 2026-07-28.
# Original location was the ephemeral session scratchpad, which was wiped.
# Session c56eada6-56b5-4b15-83cd-c6ed69cf48d9, written 07-28T08:36:13.

import numpy as np, xarray as xr, warnings, cfgrib; warnings.filterwarnings('ignore')
BOX={'Siberia':(55,75,60,180),'E. Siberia':(55,75,90,160),'Boreal':(55,70,-180,180)}
ACC=3600.0
RT='/work/bb1469/a270270/runtime/awiesm3-v3.4/'
LSMg=xr.open_dataset(RT+'Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')['lsm'].isel(time_counter=0).values
lsm_cell=xr.open_dataset('/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_forcing/input/oifs/95_4/lsm',engine='cfgrib')['lsm'].values>0.5
def bg(da,box):   # regular grid, JJA already selected
    la0,la1,lo0,lo1=box; lat=da.lat.values; lon=((da.lon.values+180)%360)-180
    yi=(lat>=la0)&(lat<=la1); xi=(lon>=lo0)&(lon<=lo1)
    v=da.values[:,yi,:][:,:,xi]; L=LSMg[yi,:][:,xi]>0.5
    w=np.cos(np.deg2rad(lat[yi]))[:,None]*np.ones(L.shape); w=np.where(L,w,0.)
    return (v*w).sum(axis=(1,2)).mean()/w.sum()
A='/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_base/outdata/oifs/'
FLUX={'ssr','str','ssrd','strd','sshf','slhf','tsr','ttr','e','sf','lsp','cp'}
VARS=['2t','skt','ssr','str','ssrd','strd','sshf','slhf','tcc','lcc','mcc','fal','sd','swvl1']
amip={b:{} for b in BOX}
for v in VARS:
    acc={b:[] for b in BOX}
    for y in range(1872,1880):
        ds=xr.open_dataset(A+f'atm_remapped_1m_{v}_1m_{y}-{y}.nc')
        da=ds[v].isel(time_counter=[5,6,7])          # JJA
        if v in FLUX: da=da/ACC
        for b,box in BOX.items(): acc[b].append(bg(da,box))
        ds.close()
    for b in BOX: amip[b][v]=float(np.mean(acc[b]))
# CRUNCEP3 JJA
ds=xr.open_dataset('/work/ab0995/a270270/input/cruncep_v7/CRUNCEP_noLPJG_1d_1901-1910_TCO95_calibrated_v3.nc')
lat=ds.lat.values; lon=((ds.lon.values+180)%360)-180
mon=ds['time_counter'].dt.month.values; jja=np.isin(mon,[6,7,8])
cru={b:{} for b in BOX}
for b,(la0,la1,lo0,lo1) in BOX.items():
    m=(lat>=la0)&(lat<=la1)&(lon>=lo0)&(lon<=lo1)&lsm_cell
    w=np.cos(np.deg2rad(lat[m]))
    for v in ['tas','rsns','rlns','pr','sfcWind','hurs']:
        cru[b][v]=float(np.average(ds[v].values[jja][:,m],axis=1,weights=w).mean())
ds.close()
print(f"{'':13s} {'Siberia':>22s} {'E. Siberia':>22s} {'Boreal':>22s}")
def row(name,f):
    print(f"{name:13s} "+" ".join(f"{f(b):>22s}" for b in ['Siberia','E. Siberia','Boreal']))
row('T2m  A/C/d',lambda b: f"{amip[b]['2t']-273.15:6.2f}/{cru[b]['tas']-273.15:6.2f}/{amip[b]['2t']-cru[b]['tas']:+5.2f}")
row('SWnet A/C/d',lambda b: f"{amip[b]['ssr']:6.1f}/{cru[b]['rsns']:6.1f}/{amip[b]['ssr']-cru[b]['rsns']:+5.1f}")
row('LWnet A/C/d',lambda b: f"{amip[b]['str']:6.1f}/{cru[b]['rlns']:6.1f}/{amip[b]['str']-cru[b]['rlns']:+5.1f}")
row('SWdn (A)',lambda b: f"{amip[b]['ssrd']:6.1f}")
row('albedo (A)',lambda b: f"{1-amip[b]['ssr']/amip[b]['ssrd']:6.3f}  fal={amip[b]['fal']:.3f}")
row('LWdn (A)',lambda b: f"{amip[b]['strd']:6.1f}")
row('SH/LH (A)',lambda b: f"{amip[b]['sshf']:6.1f}/{amip[b]['slhf']:6.1f}  B={amip[b]['sshf']/amip[b]['slhf']:.2f}")
row('cloud t/l/m',lambda b: f"{amip[b]['tcc']:.3f}/{amip[b]['lcc']:.3f}/{amip[b]['mcc']:.3f}")
row('skt-2t (A)',lambda b: f"{amip[b]['skt']-amip[b]['2t']:+6.2f}")
row('swvl1,sd(A)',lambda b: f"{amip[b]['swvl1']:.3f} {amip[b]['sd']*1000:.1f}mm")
row('pr A/C',lambda b: f"{'n/a':>6s}/{cru[b]['pr']*86400:6.2f} mm/d")
row('wind A/C',lambda b: f"{'n/a':>6s}/{cru[b]['sfcWind']:6.2f}")
