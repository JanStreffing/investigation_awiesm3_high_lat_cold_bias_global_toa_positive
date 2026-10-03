"""Sea ice in AWI-ESM2 PI against AWI-ESM3 (11X, 13A), identical definitions on each mesh.

AWI-ESM2's FESOM-2.6 output carries no coordinates (126858 nodes, CORE2), so node lat
and area come from the mesh files: area = 1/3 of the adjacent triangle areas, computed
from 3-D unit vectors so longitude wrap cannot bite.  AWI-ESM3 uses its CORE3 mesh the
same way, so extent and volume are computed identically on both.

AWI-ESM2 reference: runtime/AWIESM2_pict, linked from a colleague's
/work/ba1066/a270107/esm_tools/EXP/PI_wisofix_c (mesh from its namelist.config).

Metrics: extent (a_ice >= 0.15 area), NH March volume and floe thickness (ratio of
means), snow on ice (ratio of band means, 60-90N).  TRAPS: a_ice DAILY, m_ice/m_snow
MONTHLY; mean(m_snow/a_ice) over nodes is unreliable, use ratio of means.

Usage:  python3 scripts/analysis/awiesm2_vs_awiesm3_ice.py
"""
import os, glob
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
RE=6.371e6
_ML=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
R3='/work/bb1469/a270092/runtime/awiesm3-v3.4'
A2='/work/ab0246/a270092/runtime/AWIESM2_pict/outdata/fesom'

def meshpath_from(cfgglob):
    for c in sorted(glob.glob(cfgglob)):
        for line in open(c):
            if 'meshpath' in line.lower():
                return line.split('=')[1].strip().strip("'\", ").rstrip('/')+'/'
def read_mesh(mp):
    with open(mp+'nod2d.out') as f:
        n=int(f.readline()); xy=np.loadtxt(f,max_rows=n,usecols=(1,2))
    with open(mp+'elem2d.out') as f:
        ne=int(f.readline()); el=np.loadtxt(f,max_rows=ne,dtype=np.int64)-1
    lon,lat=np.deg2rad(xy[:,0]),np.deg2rad(xy[:,1])
    P=np.stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)],1)
    a,b,c=P[el[:,0]],P[el[:,1]],P[el[:,2]]
    ta=0.5*np.linalg.norm(np.cross(b-a,c-a),axis=1)*RE**2
    area=np.zeros(n); np.add.at(area,el.ravel(),np.repeat(ta/3.0,3))
    return xy[:,1],area

mp2=meshpath_from('/work/ba1066/a270107/esm_tools/EXP/PI_wisofix_c/config/fesom/namelist.config*')
if not mp2 or not os.path.exists(mp2+'nod2d.out'): mp2='/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2/'
mp3=meshpath_from(f'{R3}/13A/run_13500101-13591231/work/namelist.config')
print('AWI-ESM2 mesh:',mp2); print('AWI-ESM3 mesh:',mp3)
MESH={'CORE2':read_mesh(mp2),'CORE3':read_mesh(mp3)}
for k,(lat,ar) in MESH.items():
    print(f'  {k}: {lat.size} nodes, ocean area {ar.sum()/1e12:.1f} M km2')

def monthly(a):
    n=a.shape[0]
    if n==12: return a
    ml=_ML.copy(); ml[1]=29 if n==366 else 28
    e=np.cumsum(ml); s=np.r_[0,e[:-1]]
    return np.stack([np.nanmean(a[s[m]:e[m]],0) for m in range(12)])

RUNS=[('AWI-ESM2 PI','CORE2',lambda v,y:f'{A2}/{v}.fesom.{y}.nc',range(5991,6001)),
      ('11X','CORE3',lambda v,y:f'{R3}/11X/outdata/fesom/{v}.fesom.{y}.nc',range(1350,1360)),
      ('13A','CORE3',lambda v,y:f'{R3}/13A/outdata/fesom/{v}.fesom.{y}.nc',range(1350,1360))]
rows=[]
for name,mesh,pf,years in RUNS:
    lat,ar=MESH[mesh]; acc=[]
    for y in years:
        if not os.path.exists(pf('a_ice',y)): continue
        d={}
        for v in ('a_ice','m_ice','m_snow'):
            with xr.open_dataset(pf(v,y),decode_times=False) as ds:
                x=monthly(np.nan_to_num(np.asarray(ds[v].values,dtype=np.float64)))
            assert x.shape[1]==lat.size, (name,v,x.shape,lat.size)
            d[v]=x
        nh,sh,arc=lat>0,lat<0,lat>=60
        ext=lambda k,m: float(ar[m&(d['a_ice'][k]>=0.15)].sum())/1e12
        vol=lambda k,m: float((d['m_ice'][k][m]*ar[m]).sum())/1e12
        snf=lambda k: float((d['m_snow'][k][arc]*ar[arc]).sum()/max((d['a_ice'][k][arc]*ar[arc]).sum(),1))
        acc.append([ext(2,nh),ext(8,nh),ext(1,sh),ext(8,sh),vol(2,nh),
                    vol(2,nh)/max(float((d['a_ice'][2][nh]*ar[nh]).sum())/1e12,1e-9),
                    snf(4),snf(5),snf(6),snf(7)])
    rows.append((name,len(acc),np.mean(acc,0) if acc else None))

H=['NH Mar ext','NH Sep ext','SH Feb ext','SH Sep ext','NH Mar vol','NH Mar thk',
   'snow May','snow Jun','snow Jul','snow Aug']
U=['M km2']*4+['1e3 km3','m']+['m']*4
print(f"\n{'':>12} "+' '.join(f'{h:>10}' for h in H)); print(f"{'':>12} "+' '.join(f'{u:>10}' for u in U))
for name,n,r in rows:
    print(f'{name:>12} '+(' '.join(f'{x:10.2f}' for x in r) if r is not None else 'missing')+f'   ({n} yr)')
print(f"{'obs':>12} {15.0:10.2f} {'~7':>10} {'~3':>10} {18.5:10.2f}   (NH Mar/SH Sep as used in this campaign)")
