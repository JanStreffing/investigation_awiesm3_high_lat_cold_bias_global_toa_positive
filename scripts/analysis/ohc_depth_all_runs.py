"""Global ocean heat content by depth range, every 5 years from each run's start, all coupled runs.

Same method as so_heat_depth_4yr.py: annual-mean temp, rho*cp*T*cell_area*dz over the 47
FESOM levels, below-bottom points masked where temp is NaN or exactly 0.  The mesh
(core3 or core3_beta) is chosen by the node count of each file.  Runs are those in
data/toa_annual_all.csv (>= 10 complete years), minus the CMIP7 spin-ups.
Writes data/ohc_depth_all_runs.csv: run, year, J in 0-100, 100-700, 700-2000, >2000 m.

Also reads each run's last complete year.  Rows already in the CSV are reused.

Usage:  python3 scripts/analysis/ohc_depth_all_runs.py      (NPROC=8, STEP=5 by default)
"""
import os, glob, csv
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
from multiprocessing import Pool
warnings.filterwarnings('ignore')
REPO=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT={'a270092':'/work/bb1469/a270092/runtime/awiesm3-v3.4','a270270':'/work/bb1469/a270270/runtime/awiesm3-v3.4'}
MESH={220509:'/work/ab0246/a270092/input/fesom2/core3/mesh.nc',211567:'/work/ab0246/a270092/input/fesom2/core3_beta/mesh.nc'}
RHO,CP=1025.0,3992.0
DEPTHS=[(0,100),(100,700),(700,2000),(2000,1e9)]
STEP=int(os.environ.get('STEP',5))
GEOM={}
for n,p in MESH.items():
    with xr.open_dataset(p,decode_times=False) as m:
        GEOM[n]=(m['cell_area'].values.astype('f8'), np.abs(m['depth_bnds'].values.astype('f8')))

def tpath(run,root,y):
    base=f'{ROOT[root]}/{run}'
    for p in (f'{base}/outdata/fesom/temp.fesom.{y}.nc', *sorted(glob.glob(f'{base}/run_*/work/temp.fesom_{y}-{y}.nc'))):
        if os.path.isfile(p): return p
def one(job):
    run,root,y,p=job
    try:
        with xr.open_dataset(p,decode_times=False) as d:
            v=d['temp']
            if v.sizes.get('time',12)!=12: return None
            zd=[c for c in v.dims if c.startswith('nz')][0]
            a=v.mean('time').transpose(zd,...).values.astype(np.float64)
        a[a==0.0]=np.nan
        carea,db=GEOM[a.shape[1]]; nl=a.shape[0]
        dz=np.diff(db)[:nl]; zm=0.5*(db[:-1]+db[1:])[:nl]
        h=RHO*CP*np.where(np.isfinite(a),a,0.0)*carea[None,:]*dz[:,None]
        return (run,y,*[float(h[(zm>=d0)&(zm<d1)].sum()) for d0,d1 in DEPTHS])
    except Exception:
        return None

runs={}
for r in csv.DictReader(open(os.path.join(REPO,'data','toa_annual_all.csv'))):
    if r['run'].startswith('AWIESM7') or r['run'] in runs: continue
    hit=[k for k,v in ROOT.items() if os.path.isdir(f"{v}/{r['run']}")]
    if hit: runs[r['run']]=hit[0]
out=os.path.join(REPO,'data','ohc_depth_all_runs.csv')
old=list(csv.reader(open(out)))[1:] if os.path.isfile(out) else []   # cache: only new years are read
done={(r[0],int(r[1])) for r in old}
jobs=[]
for run,root in sorted(runs.items()):
    ys=sorted(int(os.path.basename(p).split('.')[-2]) for p in glob.glob(f'{ROOT[root]}/{run}/outdata/fesom/temp.fesom.[0-9]*.nc'))
    ys+= [int(os.path.basename(p).split('_')[-1].split('-')[0]) for p in glob.glob(f'{ROOT[root]}/{run}/run_*/work/temp.fesom_*-*.nc')]
    ys=sorted(set(ys))
    if not ys: continue
    want=sorted(set(range(ys[0], ys[-1]+1, STEP)) | {ys[-1]})   # every STEP years, plus the last
    for y in want:
        if (run,y) in done: continue
        p=tpath(run,root,y)
        if p: jobs.append((run,root,y,p))
print(f'{len(runs)} runs, {len(jobs)} run-years to read',flush=True)
with Pool(int(os.environ.get('NPROC',8))) as pool:
    res=[r for r in pool.imap_unordered(one,jobs,chunksize=1) if r]
res=sorted([(r[0],int(r[1]),*map(float,r[2:])) for r in old]+res)
with open(out,'w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['run','year','J_0_100','J_100_700','J_700_2000','J_gt2000']); w.writerows(res)
print(f'wrote {len(res)} rows to {out}')
