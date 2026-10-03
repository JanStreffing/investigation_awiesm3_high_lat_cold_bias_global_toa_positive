"""Annual global net TOA for every coupled AWI-ESM3 run with >= 10 complete years; cached as CSV.

Net TOA = tsr + ttr (accumulated J/m2 per hourly step, /3600), cos-lat weighted global mean
(weights from each file's own lat), annual mean over 12 complete months.  Looks in
outdata/oifs first, then run_*/work.  Writes data/toa_annual_all.csv (run,root,year,toa).

Usage:  python3 scripts/analysis/toa_annual_all_runs.py      (NPROC=8 by default)
"""
import os, glob, csv, sys
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
from multiprocessing import Pool
warnings.filterwarnings('ignore')
ROOTS=['/work/bb1469/a270092/runtime/awiesm3-v3.4','/work/bb1469/a270270/runtime/awiesm3-v3.4']
REPO=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SKIP=('failed_','HP_','Test_Historical')
ACC=3600.0

def files(d):
    out={}
    for p in glob.glob(f'{d}/outdata/oifs/atm_remapped_1m_tsr_*.nc')+glob.glob(f'{d}/run_*/work/atm_remapped_1m_tsr_*.nc'):
        y=int(os.path.basename(p).split('_')[-1].split('-')[0]); q=p.replace('_tsr_','_ttr_')
        if os.path.isfile(q) and y not in out: out[y]=(p,q)
    return out
def ld(p):
    with xr.open_dataset(p,decode_times=False) as d:
        k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values), np.squeeze(d['lat'].values)
def one(job):
    run,root,y,p,q=job
    try:
        s,lat=ld(p); t,_=ld(q)
        if s.shape[0]!=12 or t.shape[0]!=12: return None
        f=(s+t).mean(0)/ACC; w=np.cos(np.deg2rad(np.broadcast_to(lat[:,None],f.shape)))
        return (run,root,y,float((f*w).sum()/w.sum()))
    except Exception as e:
        return None
jobs=[]
for root in ROOTS:
    for d in sorted(glob.glob(f'{root}/*/')):
        run=os.path.basename(d.rstrip('/'))
        if run.startswith(SKIP) or any(x in run for x in ('.aborted','.failed','.cancelled','_dead_','.bak')): continue
        fs=files(d.rstrip('/'))
        if len(fs)<10: continue
        jobs+= [(run,root.split('/')[3],y,*fs[y]) for y in sorted(fs)]
print(f'{len(jobs)} run-years to read',flush=True)
with Pool(int(os.environ.get('NPROC',8))) as pool:
    res=[r for r in pool.imap_unordered(one,jobs,chunksize=4) if r]
res.sort()
out=os.path.join(REPO,'data','toa_annual_all.csv')
with open(out,'w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['run','root','year','toa']); w.writerows(res)
print(f'wrote {len(res)} rows to {out}')
