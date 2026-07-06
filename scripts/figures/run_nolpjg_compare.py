#!/usr/bin/env python
# Decisive test (a la EC-Earth): coupled WITH LPJG (Tuning_06 Baseline, core3) vs
# WITHOUT LPJG (AWI-CM3 TUNE42PI_FES27, HTESSEL + prescribed satellite veg, core2),
# same TCO95 atmosphere, same years 1370-1379. Compare boreal 2m T, surface albedo
# (1-SSR/SSRD, incl. spring) and high-veg cover cvh.
import numpy as np, xarray as xr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
T='/work/bb1469/a270092/eval/nolpjg'; EVAL='/work/bb1469/a270092/eval'
def op(tag,v):
    ds=xr.open_dataset(f'{T}/{tag}_mon_{v}.nc'); return ds[v]
def cvh(tag):
    ds=xr.open_dataset(f'{T}/{tag}_cvh.nc'); return ds['cvh'].squeeze()
lat=op('WITH','2t').lat; lon=op('WITH','2t').lon
wlat=np.cos(np.deg2rad(lat))
def boxmean(da):  # Siberia + N.America boreal land, 55-70N
    sib=da.where((lat>=55)&(lat<=70)&(lon>=60)&(lon<=140))
    can=da.where((lat>=55)&(lat<=70)&(lon>=235)&(lon<=300))
    return 0.5*(float(sib.weighted(wlat).mean(('lat','lon')))+float(can.weighted(wlat).mean(('lat','lon'))))
def alb(tag):  # monthly surface albedo
    return (1-op(tag,'ssr')/op(tag,'ssrd')).where(op(tag,'ssrd')>5)
mon=np.arange(1,13)
# seasonal cycles
T_with=[boxmean(op('WITH','2t').isel(time_counter=m) if 'time_counter' in op('WITH','2t').dims else op('WITH','2t').isel(time=m)) for m in range(12)]
def monthly_box(tag,field):
    da=field(tag); td=[d for d in da.dims if 'time' in d][0]
    return [boxmean(da.isel({td:m})) for m in range(12)]
Tw=monthly_box('WITH',lambda t:op(t,'2t')); Two=monthly_box('WITHOUT',lambda t:op(t,'2t'))
Aw=monthly_box('WITH',alb); Awo=monthly_box('WITHOUT',alb)
Tw=np.array(Tw)-273.15; Two=np.array(Two)-273.15; Aw=np.array(Aw); Awo=np.array(Awo)

# ---- OIFS-AMIP boreal 2m-T seasonal cycle (observed-SST, no-LPJG; the forcing
# generator for the offline re-spin). tas is a STATE field (unaffected by the
# flux-accumulation bug that only hit rsns/rlns/pr), so it is used as-is. No
# ssrd is emitted, so AMIP appears on the T panel only, not the albedo panel. ----
import cftime as _cf
_am=xr.open_dataset('/work/ab0246/a270092/input/lpj-guess/oifs_forcing/'
                    'AMIP_noLPJG_1d_1870-1879_TCO95_PI.nc',decode_times=False)
_alat=_am['lat'].values; _alon=_am['lon'].values
_tc=_am['time_counter']
_mon_of=np.array([d.month for d in _cf.num2date(_tc.values,_tc.attrs['units'],_tc.attrs.get('calendar','gregorian'))])
_box=(((_alat>=55)&(_alat<=70))&(((_alon>=60)&(_alon<=140))|((_alon>=235)&(_alon<=300))))
_idx=np.where(_box)[0]; _w=np.cos(np.deg2rad(_alat))[_idx]
_tb=_am['tas'].values[:,_idx]                                    # (3652, nbox)
Tam=np.array([np.average(_tb[_mon_of==m].mean(axis=0),weights=_w) for m in range(1,13)])-273.15

# ---- Fig 1: seasonal cycle ----
fig,ax=plt.subplots(1,2,figsize=(13,5))
ax[0].plot(mon,Two,'-o',color='#2e7d32',label='WITHOUT LPJG (HTESSEL, coupled SST)')
ax[0].plot(mon,Tw,'-o',color='#c62828',label='WITH LPJG (coupled)')
ax[0].plot(mon,Tam,'-s',color='#1565c0',label='OIFS-AMIP (obs SST, no LPJG)')
ax[0].set_title('Boreal 2 m temperature (55–70°N, Siberia+N.America)',fontweight='bold')
ax[0].set_xlabel('month'); ax[0].set_ylabel('°C'); ax[0].grid(alpha=.3); ax[0].legend()
ax2=ax[0].twinx(); ax2.bar(mon,Tw-Two,alpha=.2,color='b'); ax2.set_ylabel('WITH−WITHOUT [K]',color='b')
ax[1].plot(mon,Awo,'-o',color='#2e7d32',label='WITHOUT LPJG')
ax[1].plot(mon,Aw,'-o',color='#c62828',label='WITH LPJG')
ax[1].set_title('Boreal surface albedo (1−SSR/SSRD)',fontweight='bold')
ax[1].set_xlabel('month'); ax[1].set_ylabel('albedo'); ax[1].grid(alpha=.3); ax[1].legend()
plt.tight_layout(); plt.savefig(f'{EVAL}/plots/nolpjg_seasonal_cycle.png',dpi=150,bbox_inches='tight')

# ---- Fig 2: NH maps annual T, spring albedo, cvh ----
def ann(tag,v):
    da=op(tag,v); td=[d for d in da.dims if 'time' in d][0]; return da.mean(td)
def seas(tag,vfun,months):
    da=vfun(tag); td=[d for d in da.dims if 'time' in d][0]; return da.isel({td:[m-1 for m in months]}).mean(td)
rows=[('2 m T annual [°C]', ann('WITH','2t')-273.15, ann('WITHOUT','2t')-273.15, (-30,5,'turbo'),(-8,8,'RdBu_r')),
      ('Surface albedo, spring MAM', seas('WITH',alb,[3,4,5]), seas('WITHOUT',alb,[3,4,5]), (0.2,0.8,'viridis'),(-0.25,0.25,'RdBu_r')),
      ('High-veg cover cvh', cvh('WITH'), cvh('WITHOUT'), (0,0.8,'YlGn'),(-0.5,0.5,'BrBG'))]
fig=plt.figure(figsize=(15,12)); proj=ccrs.NorthPolarStereo()
for r,(name,w,wo,(vmn,vmx,cm),(dmn,dmx,dcm)) in enumerate(rows):
    for c,(t,z,vmn2,vmx2,cm2) in enumerate([(f'{name}\nWITH LPJG',w,vmn,vmx,cm),
                                            (f'{name}\nWITHOUT LPJG',wo,vmn,vmx,cm),
                                            (f'{name}\nWITH − WITHOUT',w-wo,dmn,dmx,dcm)]):
        ax=fig.add_subplot(3,3,r*3+c+1,projection=proj); ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
        im=ax.pcolormesh(lon,lat,np.ma.masked_invalid(z.values),cmap=cm2,vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),shading='auto')
        ax.add_feature(cfeature.OCEAN,facecolor='lightblue',alpha=.5,zorder=2); ax.coastlines(linewidth=.4,zorder=3)
        ax.set_title(t,fontsize=9.5,fontweight='bold'); plt.colorbar(im,ax=ax,shrink=.6,pad=.03)
fig.suptitle('AWI-ESM3 WITH LPJG vs AWI-CM3 WITHOUT LPJG (HTESSEL), TCO95, yr1370-79',fontsize=13,fontweight='bold',y=1.0)
plt.tight_layout(); plt.savefig(f'{EVAL}/plots/nolpjg_maps.png',dpi=145,bbox_inches='tight')

print(f"=== boreal box (55-70N, Siberia+Canada) ===")
print(f"annual 2m T   WITH={np.mean(Tw):.2f}  WITHOUT={np.mean(Two):.2f}  AMIP={np.mean(Tam):.2f}  diff(W-WO)={np.mean(Tw)-np.mean(Two):+.2f} C")
print(f"DJF    2m T   WITH={np.mean(Tw[[11,0,1]]):.2f}  WITHOUT={np.mean(Two[[11,0,1]]):.2f}  AMIP={np.mean(Tam[[11,0,1]]):.2f}  diff(W-WO)={np.mean(Tw[[11,0,1]])-np.mean(Two[[11,0,1]]):+.2f} C")
print(f"MAM    2m T   WITH={np.mean(Tw[2:5]):.2f}  WITHOUT={np.mean(Two[2:5]):.2f}  AMIP={np.mean(Tam[2:5]):.2f}  diff(W-WO)={np.mean(Tw[2:5])-np.mean(Two[2:5]):+.2f} C")
print(f"WITH-WITHOUT by month [K]: "+" ".join(f"{m+1}:{Tw[m]-Two[m]:+.1f}" for m in range(12)))
print(f"MAM albedo    WITH={np.nanmean(Aw[2:5]):.3f}  WITHOUT={np.nanmean(Awo[2:5]):.3f}  diff={np.nanmean(Aw[2:5])-np.nanmean(Awo[2:5]):+.3f}")
print(f"JJA albedo    WITH={np.nanmean(Aw[5:8]):.3f}  WITHOUT={np.nanmean(Awo[5:8]):.3f}  diff={np.nanmean(Aw[5:8])-np.nanmean(Awo[5:8]):+.3f}")
print(f"cvh           WITH={boxmean(cvh('WITH')):.3f}  WITHOUT={boxmean(cvh('WITHOUT')):.3f}  diff={boxmean(cvh('WITH'))-boxmean(cvh('WITHOUT')):+.3f}")
print("saved nolpjg_seasonal_cycle.png, nolpjg_maps.png")
