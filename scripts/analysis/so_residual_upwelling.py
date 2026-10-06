"""Resolved, GM-bolus and residual vertical velocity in the newly convecting Southern Ocean region
of PICAL_crunveg_kpp_gm1000 (region NEW of so_destabilisation.py), arm and TKE-only line.
Usage: so_residual_upwelling.py   (reval environment)"""
import numpy as np, xarray as xr
P='/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
R={'arm':f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_kpp_gm1000/outdata/fesom','tke':f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082/outdata/fesom'}
md=xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc'); lat=md['lat'].values; zi=np.abs(md['nz'].values); area=md['nod_area'].values[0]
od=(md['ulevels_nod2D'].values==1)&(np.abs(md['zbar_n_bottom'].values)>2000)&(lat<-55)
m=lambda r:np.abs(xr.open_dataset(f'{R[r]}/MLD2.fesom.2177.nc')['MLD2'].isel(time=8).values)
reg=od&(m('arm')>1000)&(m('tke')<500); idx=np.where(reg)[0]; w=area[idx]/area[idx].sum(); Y=86400*365
print('region NEW, annual mean vertical velocity [m/yr], positive up: resolved w, GM bolus w, sum')
for z in (100,200,300):
    k=int(np.argmin(np.abs(zi-z)))
    for r in R:
        for y in (2170,2171,2173):
            a=float(np.nansum(xr.open_dataset(f'{R[r]}/w.fesom.{y}.nc')['w'].isel(nz1=k).mean('time').values[idx]*w))*Y
            b=float(np.nansum(xr.open_dataset(f'{R[r]}/bolus_w.fesom.{y}.nc')['bolus_w'].isel(nz1=k).mean('time').values[idx]*w))*Y
            print(f'{z:4d} m {r} {y}: w {a:6.1f}  bolus {b:6.1f}  residual {a+b:6.1f}')
