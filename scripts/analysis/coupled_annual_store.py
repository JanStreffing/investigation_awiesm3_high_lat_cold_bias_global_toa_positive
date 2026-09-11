"""Annual diagnostics for every coupled AWI-ESM3 run, kept in one NetCDF store that only grows.

data/coupled_annual_diag.nc, dims (run, year), year = model calendar year:
  toa          global net TOA [W m-2]: tsr + ttr (accumulated J/m2 per hourly step, /3600),
               cos-lat weighted global mean over 12 complete months
  ohc_0_100, ohc_100_700, ohc_700_2000, ohc_gt2000
               global ocean heat content [J] in depth ranges: annual-mean temp,
               rho*cp*T*cell_area*dz over the 47 FESOM levels, below-bottom points masked where
               temp is exactly 0 (as so_heat_depth_4yr.py); mesh chosen by node count
  root(run)    runtime root the run was found in;  nodes(run) FESOM mesh size (0 = unknown)

Run it again to add new runs or years: only missing (run, year, quantity) entries are computed,
and the store is rewritten every CHECKPOINT results, so an interrupted run loses little.
On the first run it seeds from data/toa_annual_all.csv and data/ohc_depth_all_runs.csv.

Usage:  python3 scripts/analysis/coupled_annual_store.py
        RUNS=15F,16A ONLY=ohc NPROC=8 python3 scripts/analysis/coupled_annual_store.py
"""
import os, glob, csv
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
from multiprocessing import Pool
warnings.filterwarnings('ignore')

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STORE = os.path.join(REPO, 'data', 'coupled_annual_diag.nc')
ROOTS = {'a270092': '/work/bb1469/a270092/runtime/awiesm3-v3.4',
         'a270270': '/work/bb1469/a270270/runtime/awiesm3-v3.4'}
SKIP_PREFIX = ('failed_', 'HP_', 'Test_Historical', 'AWIESM7')
SKIP_PART = ('.aborted', '.failed', '.cancelled', '_dead_', '.bak')
MESH = {220509: '/work/ab0246/a270092/input/fesom2/core3/mesh.nc',
        211567: '/work/ab0246/a270092/input/fesom2/core3_beta/mesh.nc'}
RHO, CP, ACC = 1025.0, 3992.0, 3600.0
DEPTHS = [('ohc_0_100', 0, 100), ('ohc_100_700', 100, 700), ('ohc_700_2000', 700, 2000),
          ('ohc_gt2000', 2000, 1e9)]
OHC = [d[0] for d in DEPTHS]
VARS = ['toa'] + OHC
CHECKPOINT = int(os.environ.get('CHECKPOINT', 200))
GEOM = {}


def geom(n):
    if n not in GEOM:
        with xr.open_dataset(MESH[n], decode_times=False) as m:
            GEOM[n] = (m['cell_area'].values.astype('f8'), np.abs(m['depth_bnds'].values.astype('f8')))
    return GEOM[n]


def files(base, pattern_out, pattern_work, year_of):
    out = {}
    for p in sorted(glob.glob(f'{base}/outdata/{pattern_out}')) + sorted(glob.glob(f'{base}/run_*/work/{pattern_work}')):
        try:
            y = year_of(os.path.basename(p))
        except ValueError:
            continue
        out.setdefault(y, p)
    return out


def toa_files(base):
    f = files(base, 'oifs/atm_remapped_1m_tsr_*.nc', 'atm_remapped_1m_tsr_*.nc',
              lambda b: int(b.split('_')[-1].split('-')[0]))
    return {y: (p, p.replace('_tsr_', '_ttr_')) for y, p in f.items() if os.path.isfile(p.replace('_tsr_', '_ttr_'))}


def temp_files(base):
    return files(base, 'fesom/temp.fesom.[0-9]*.nc', 'temp.fesom_*-*.nc',
                 lambda b: int(b.split('.')[-2]) if b.startswith('temp.fesom.') else int(b.split('_')[-1].split('-')[0]))


def load(p):
    with xr.open_dataset(p, decode_times=False) as d:
        k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values), np.squeeze(d['lat'].values)


def work(job):
    kind, run, y, p = job
    try:
        if kind == 'toa':
            (s, lat), (t, _) = load(p[0]), load(p[1])
            if s.shape[0] != 12 or t.shape[0] != 12:
                return None
            f = (s + t).mean(0) / ACC
            w = np.cos(np.deg2rad(np.broadcast_to(lat[:, None], f.shape)))
            return (run, y, {'toa': float((f * w).sum() / w.sum())}, 0)
        with xr.open_dataset(p, decode_times=False) as d:
            v = d['temp']
            if v.sizes.get('time', 12) != 12:
                return None
            zd = [c for c in v.dims if c.startswith('nz')][0]
            a = v.mean('time').transpose(zd, ...).values.astype(np.float64)
        a[a == 0.0] = np.nan
        carea, db = geom(a.shape[1]); nl = a.shape[0]
        dz = np.diff(db)[:nl]; zm = 0.5 * (db[:-1] + db[1:])[:nl]
        h = RHO * CP * np.where(np.isfinite(a), a, 0.0) * carea[None, :] * dz[:, None]
        return (run, y, {n: float(h[(zm >= d0) & (zm < d1)].sum()) for n, d0, d1 in DEPTHS}, a.shape[1])
    except Exception:
        return None


# ---- existing store, or seed from the CSV caches -------------------------------------------
VAL = {v: {} for v in VARS}          # VAL[var][(run, year)] = value
META = {}                            # META[run] = [root, nodes]
if os.path.isfile(STORE):
    with xr.open_dataset(STORE) as ds:
        for i, run in enumerate(ds['run'].values.astype(str)):
            META[run] = [str(ds['root'].values[i]), int(ds['nodes'].values[i])]
            for v in VARS:
                col = ds[v].values[i]
                for j, y in enumerate(ds['year'].values):
                    if np.isfinite(col[j]):
                        VAL[v][(run, int(y))] = float(col[j])
else:
    t = os.path.join(REPO, 'data', 'toa_annual_all.csv')
    if os.path.isfile(t):
        for r in csv.DictReader(open(t)):
            VAL['toa'][(r['run'], int(r['year']))] = float(r['toa'])
    o = os.path.join(REPO, 'data', 'ohc_depth_all_runs.csv')
    if os.path.isfile(o):
        for r in csv.DictReader(open(o)):
            for v, c in zip(OHC, ('J_0_100', 'J_100_700', 'J_700_2000', 'J_gt2000')):
                VAL[v][(r['run'], int(r['year']))] = float(r[c])


def write():
    runs = sorted({str(k[0]) for v in VARS for k in VAL[v]} | {str(r) for r in META})
    years = sorted({k[1] for v in VARS for k in VAL[v]})
    ri = {r: i for i, r in enumerate(runs)}; yi = {y: j for j, y in enumerate(years)}
    data = {}
    for v in VARS:
        a = np.full((len(runs), len(years)), np.nan)
        for (r, y), x in VAL[v].items():
            a[ri[str(r)], yi[int(y)]] = x
        data[v] = (('run', 'year'), a)
    data['root'] = (('run',), np.array([META.get(r, ['?', 0])[0] for r in runs], dtype=object))
    data['nodes'] = (('run',), np.array([META.get(r, ['?', 0])[1] for r in runs], dtype=np.int32))
    ds = xr.Dataset(data, coords={'run': np.array(runs, dtype=str), 'year': np.array(years, dtype=np.int32)})
    ds['toa'].attrs = {'units': 'W m-2', 'long_name': 'global net TOA, tsr+ttr, annual mean'}
    for v in OHC:
        ds[v].attrs = {'units': 'J', 'long_name': f'global ocean heat content, {v[4:].replace("_", "-")} m'}
    ds.attrs = {'description': 'annual diagnostics of the AWI-ESM3 coupled runs; grown by '
                               'scripts/analysis/coupled_annual_store.py'}
    tmp = STORE + '.tmp'
    ds.to_netcdf(tmp); os.replace(tmp, STORE)
    return len(runs), len(years)


# ---- find what is missing ----------------------------------------------------------------
only = os.environ.get('ONLY', '')
want_runs = set(filter(None, os.environ.get('RUNS', '').split(',')))
jobs = []
for root, base_root in ROOTS.items():
    for d in sorted(glob.glob(f'{base_root}/*/')):
        run = os.path.basename(d.rstrip('/'))
        if run.startswith(SKIP_PREFIX) or any(x in run for x in SKIP_PART):
            continue
        if want_runs and run not in want_runs:
            continue
        base = d.rstrip('/')
        tf = toa_files(base)
        if len(tf) < 5 and run not in want_runs:      # not a coupled production run
            continue
        META.setdefault(run, [root, 0])[0] = root
        if only in ('', 'toa'):
            jobs += [('toa', run, y, p) for y, p in tf.items() if (run, y) not in VAL['toa']]
        if only in ('', 'ohc'):
            jobs += [('ohc', run, y, p) for y, p in temp_files(base).items() if (run, y) not in VAL['ohc_700_2000']]
print(f'{len(jobs)} missing entries ({sum(j[0]=="toa" for j in jobs)} toa, {sum(j[0]=="ohc" for j in jobs)} ohc)', flush=True)

n = 0
if jobs:
    with Pool(int(os.environ.get('NPROC', 8))) as pool:
        for res in pool.imap_unordered(work, jobs, chunksize=1):
            if res is None:
                continue
            run, y, vals, nodes = res
            for v, x in vals.items():
                VAL[v][(run, y)] = x
            if nodes:
                META.setdefault(run, ['?', 0])[1] = nodes
            n += 1
            if n % CHECKPOINT == 0:
                print(f'  checkpoint: {n} entries, store {write()}', flush=True)
nr, ny = write()
print(f'done: {n} new entries; store {STORE} holds {nr} runs x {ny} years')
