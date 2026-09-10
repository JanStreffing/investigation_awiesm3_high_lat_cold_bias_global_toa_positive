"""Detection thresholds for 4-year coupled means, from 11X's 40-year interannual scatter.

For each metric: per-year values 1350-1389, linear trend removed (the run drifts), sd of
the residuals, and the 95 % threshold for a difference of two independent 4-year means,
1.96 * sd * sqrt(2/4).  Differences between arms below it are not resolved.
Same definitions as coupled_state_vs_obs.py (OIFS remapped, cos-lat weights, ERA5-free).

Usage:  python3 scripts/analysis/detection_threshold_4yr.py
        ARM=15F Y1=1359 SERIES=1 python3 ...   (shorter run; also print the per-year values)
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'; ARM=os.environ.get('ARM','11X'); ACC=3600.0
Y0,Y1=1350,int(os.environ.get('Y1',1389)); SERIES=os.environ.get('SERIES','0')=='1'
LSMF=('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
      'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
with xr.open_dataset(LSMF,decode_times=False) as d:
    m=np.squeeze(d['lsm'].values); m=m[0] if m.ndim==3 else m
    lat=np.squeeze(d['lat'].values); lon=np.squeeze(d['lon'].values)
land,ocean=m>0.5,m<=0.5; L=np.broadcast_to(lat[:,None],m.shape); W=np.cos(np.deg2rad(L))
bnd=lambda lo,hi:(L>=lo)&(L<hi); G=np.ones_like(m,bool)
sib=land&bnd(50,70)&np.broadcast_to(((lon>=60)&(lon<=140))[None,:],m.shape)
cell=W*6.371e6**2*np.deg2rad(abs(lat[1]-lat[0]))*2*np.pi/m.shape[1]
am=lambda f,s: float(np.average(f[s&np.isfinite(f)],weights=W[s&np.isfinite(f)]))
def ld(v,y):
    with xr.open_dataset(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc',decode_times=False) as d:
        k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]; return np.squeeze(d[k].values)
rows=[]
for y in range(Y0,Y1+1):
    t=ld('2t',y); tsr=ld('tsr',y); ttr=ld('ttr',y); ci=ld('ci',y)
    rows.append([am(t[[0,1,11]].mean(0),land&bnd(60,90)), am(t[[5,6,7]].mean(0),sib),
                 am(t.mean(0),G), am(t.mean(0),bnd(-90,-60)), am((tsr+ttr).mean(0)/ACC,G),
                 am((tsr-ld('tsrc',y)).mean(0)/ACC,bnd(-65,-45)),
                 float(((ci[2]>0.15)*cell)[bnd(0,90)&ocean].sum())/1e12,
                 float(((ci[8]>0.15)*cell)[bnd(0,90)&ocean].sum())/1e12,
                 float(((ci[8]>0.15)*cell)[bnd(-90,0)&ocean].sum())/1e12])
A=np.array(rows); x=np.arange(len(A))
names=['T2m DJF 60-90N land','T2m JJA Siberia','T2m ANN global','T2m ANN 90-60S','net TOA global',
       'SW CRE 45-65S','NH ice EXTENT Mar','NH ice EXTENT Sep','SH ice EXTENT Sep']
print(f'{ARM} {Y0}-{Y1}, detrended interannual sd and 95 % threshold for a 4-yr-mean difference')
for j,n in enumerate(names):
    r=A[:,j]-np.polyval(np.polyfit(x,A[:,j],1),x); sd=r.std(ddof=2)
    print(f'  {n:<22} sd {sd:7.3f}   threshold {1.96*sd*np.sqrt(2/4):7.3f}')
if SERIES:
    short=['DJF60-90N','JJASib','T2m glob','T2m 90-60S','netTOA','SWCRE SO','NHext Mar','NHext Sep','SHext Sep']
    print('\n  year '+''.join(f'{n:>11}' for n in short))
    for y,r in zip(range(Y0,Y1+1),A):
        print(f'  {y} '+''.join(f'{v:11.3f}' for v in r))
    print('  trend/decade '+''.join(f'{10*np.polyfit(x,A[:,j],1)[0]:11.3f}' for j in range(A.shape[1]))[3:])
