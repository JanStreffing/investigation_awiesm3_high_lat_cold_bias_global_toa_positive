"""How much more sea ice should an 1850 model have than the satellite-era observations?

HadISST2.2 sea-ice concentration, 1850-2019 monthly on 1x1.  Extent = area of cells with
sic >= 0.15, true spherical cell areas.  Compares 1850-1879 (the PI target window) with
1979-2008 (the satellite era the campaign's 15.0 / 18.5 M km2 references come from).

CAVEAT THAT MUST TRAVEL WITH THE NUMBER.  Before the satellite era (pre-1979) HadISST2
Arctic ice comes from sparse ship/aircraft records and reconstruction, and before ~1950
much of it is climatology-infilled.  The 1850-79 values are a reconstruction, not a
measurement; treat the PI-minus-satellite offset as an estimate with an uncertainty of
several tenths of a M km2 at least, larger in the Southern Hemisphere.

Usage:  python3 scripts/analysis/hadisst_pi_vs_satellite_extent.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
F='/work/ab0246/a270092/obs/hadisst2/HadISST.2.2.0.0_sea_ice_concentration.nc'
RE=6.371e6
with xr.open_dataset(F,decode_times=False) as d:
    sic=np.asarray(d['sic'].values,dtype=np.float64)
    lat=np.asarray(d['latitude'].values); lon=np.asarray(d['longitude'].values)
    t=np.asarray(d['time'].values); tu=d['time'].attrs.get('units')
print('time units:',tu,' n months:',len(t))
sic[~np.isfinite(sic)]=0.0; sic[sic<0]=0.0
if np.nanmax(sic)>1.5: sic=sic/100.0
dl=np.deg2rad(abs(lon[1]-lon[0])); h=abs(lat[1]-lat[0])/2
area=RE**2*dl*np.abs(np.sin(np.deg2rad(lat+h))-np.sin(np.deg2rad(lat-h)))
A=np.broadcast_to(area[:,None],sic.shape[1:])
# time axis is monthly from 1850-01 (days since 1850-1-1); index -> (year, month)
yr=1850+np.arange(len(t))//12; mo=np.arange(len(t))%12
nh=np.broadcast_to((lat>0)[:,None],A.shape); sh=~nh
def ext(i,hemi): return float(A[hemi&(sic[i]>=0.15)].sum())/1e12
def clim(y0,y1,m,hemi):
    idx=np.where((yr>=y0)&(yr<=y1)&(mo==m))[0]
    return np.mean([ext(i,hemi) for i in idx]), len(idx)
print(f"\n{'':>22} {'NH Mar':>8} {'NH Sep':>8} {'SH Feb':>8} {'SH Sep':>8}   M km2")
res={}
for lab,(y0,y1) in (('1850-1879 (PI window)',(1850,1879)),('1979-2008 (satellite)',(1979,2008))):
    r=[clim(y0,y1,2,nh)[0],clim(y0,y1,8,nh)[0],clim(y0,y1,1,sh)[0],clim(y0,y1,8,sh)[0]]
    res[lab]=r; print(f'{lab:>22} '+' '.join(f'{x:8.2f}' for x in r))
d=np.array(res['1850-1879 (PI window)'])-np.array(res['1979-2008 (satellite)'])
print(f"{'PI minus satellite':>22} "+' '.join(f'{x:+8.2f}' for x in d))
