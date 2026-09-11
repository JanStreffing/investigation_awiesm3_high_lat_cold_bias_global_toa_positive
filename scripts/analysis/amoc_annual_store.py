"""Annual AMOC indices for the coupled runs, kept in one NetCDF store that only grows.

data/amoc_annual_diag.nc, dims (run, year):
  amoc26      Atlantic overturning maximum at 26.5N below 500 m, Eulerian (from w) [Sv]
  amoc26_res  the same for the residual flow, w + bolus_w (adds the GM eddy-induced part)
  amoc4060    maximum over 40-60N below 500 m, Eulerian
  amoc4060_res  the same, residual
Method: pyfesom2 xmoc_data (as release_evaluation_tool2 part17_moc.py) on the annual-mean
field, Atlantic_MOC basin mask, 181 latitude bins, element areas and levels from the mesh's
fesom.mesh.diag.nc (core3 or core3_beta, chosen by node count).  Only missing (run, year)
entries are computed.  Needs pyfesom2, so run it in the reval environment:

  source ~/loadconda.sh && conda activate reval
  python3 scripts/analysis/amoc_annual_store.py            (RUNS=15F,16C ... NPROC=6)
"""
import os, glob
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import warnings; warnings.filterwarnings('ignore')
import numpy as np, xarray as xr
import pyfesom2 as pf
from pyfesom2.ut import get_mask
from multiprocessing import Pool

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STORE = os.path.join(REPO, 'data', 'amoc_annual_diag.nc')
ROOT = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESHDIR = {220509: '/work/ab0246/a270092/input/fesom2/core3/',
           211567: '/work/ab0246/a270092/input/fesom2/core3_beta/'}
DIRS = {'11E': 'Tuning_test_11E_swemin15_K1', '11G': 'Tuning_test_11G_inppmin50k',
        '11I': 'Tuning_test_11I_v2soil'}
RUNS = [r for r in os.environ.get('RUNS', '11X,11W,11E,11G,11I,15F,16A,16B,16C,16D').split(',') if r]
VARS = ['amoc26', 'amoc26_res', 'amoc4060', 'amoc4060_res']
NLATS = 181
G = {}


def geometry(n):
    if n not in G:
        mesh = pf.load_mesh(MESHDIR[n])
        md = xr.open_dataset(os.path.join(MESHDIR[n], 'fesom.mesh.diag.nc'))
        fx, fy = pf.diagnostics.compute_face_coords(mesh)
        G[n] = dict(mesh=mesh, el_area=md['elem_area'], nlevels=md['nlevels'].values - 1,
                    fx=fx, fy=fy, mask=get_mask(mesh, 'Atlantic_MOC'), z=np.abs(np.asarray(mesh.zlev)))
    return G[n]


def annual(p, name):
    with xr.open_dataset(p, decode_times=False) as d:
        k = name if name in d else [c for c in d.data_vars if d[c].ndim == 3][0]
        v = d[k]
        if v.sizes.get('time', 12) != 12:
            return None
        zd = [c for c in v.dims if c.startswith('nz')][0]
        return v.mean('time').transpose('nod2', zd).values.astype(np.float64)


def index(g, w, lat0, lat1):
    lats, moc = pf.xmoc_data(g['mesh'], w, nlats=NLATS, mask=g['mask'], el_area=g['el_area'],
                             nlevels=g['nlevels'], face_x=g['fx'], face_y=g['fy'])
    moc = np.ma.filled(moc, np.nan)                     # (nlats, nlev)
    zsel = g['z'][:moc.shape[1]] >= 500.0
    lsel = (lats >= lat0) & (lats <= lat1)
    if lat0 == lat1:
        lsel = np.abs(lats - lat0) == np.min(np.abs(lats - lat0))
    return float(np.nanmax(moc[np.ix_(lsel, zsel)]))


def work(job):
    run, y, pw, pb = job
    try:
        w = annual(pw, 'w')
        if w is None:
            return None
        g = geometry(w.shape[0])
        out = {'amoc26': index(g, w, 26.5, 26.5), 'amoc4060': index(g, w, 40.0, 60.0)}
        if pb:
            b = annual(pb, 'bolus_w')
            if b is not None and b.shape == w.shape:
                out['amoc26_res'] = index(g, w + b, 26.5, 26.5)
                out['amoc4060_res'] = index(g, w + b, 40.0, 60.0)
        return run, y, out
    except Exception as e:
        print(f'  {run} {y}: {e}', flush=True)
        return None


VAL = {v: {} for v in VARS}
if os.path.isfile(STORE):
    with xr.open_dataset(STORE) as ds:
        for i, r in enumerate(str(x) for x in ds['run'].values):
            for v in VARS:
                for j, y in enumerate(ds['year'].values):
                    x = ds[v].values[i, j]
                    if np.isfinite(x):
                        VAL[v][(r, int(y))] = float(x)


def write():
    runs = sorted({str(k[0]) for v in VARS for k in VAL[v]}); years = sorted({int(k[1]) for v in VARS for k in VAL[v]})
    ri = {r: i for i, r in enumerate(runs)}; yi = {y: j for j, y in enumerate(years)}
    data = {}
    for v in VARS:
        a = np.full((len(runs), len(years)), np.nan)
        for (r, y), x in VAL[v].items():
            a[ri[str(r)], yi[int(y)]] = x
        data[v] = (('run', 'year'), a, {'units': 'Sv'})
    ds = xr.Dataset(data, coords={'run': np.array(runs, dtype=str), 'year': np.array(years, dtype=np.int32)})
    ds.attrs = {'description': 'annual AMOC indices, grown by scripts/analysis/amoc_annual_store.py'}
    tmp = STORE + '.tmp'; ds.to_netcdf(tmp); os.replace(tmp, STORE)
    return len(runs), len(years)


jobs = []
for run in RUNS:
    base = f'{ROOT}/{DIRS.get(run, run)}/outdata/fesom'
    for pw in sorted(glob.glob(f'{base}/w.fesom.[0-9]*.nc')):
        y = int(os.path.basename(pw).split('.')[-2])
        if (run, y) in VAL['amoc26']:
            continue
        pb = f'{base}/bolus_w.fesom.{y}.nc'
        jobs.append((run, y, pw, pb if os.path.isfile(pb) else None))
print(f'{len(jobs)} missing run-years', flush=True)
n = 0
if jobs:
    with Pool(int(os.environ.get('NPROC', 6))) as pool:
        for res in pool.imap_unordered(work, jobs, chunksize=1):
            if res is None:
                continue
            run, y, out = res
            for v, x in out.items():
                VAL[v][(run, y)] = x
            n += 1
            if n % 50 == 0:
                print(f'  checkpoint {n}: store {write()}', flush=True)
print(f'done: {n} new; store {STORE} holds {write()}')
