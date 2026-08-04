"""Shared, cached per-run extraction of Siberian-box monthly means.

Why this exists: `eval_round10_A.py` was made fast with a per-run JSON cache
(2m09s -> 7.7s), and then three UNCACHED scripts were added to `evaluate.sh`
(`noise_floor.py`, `seasonal_by_run.py`, `monthly_lever_check.py`), each of which
re-reads 2t for every run and every year. They all want the same numbers, so the
wrapper ended up slower than running the old scripts by hand.

One cache, shared. `monthly(run, var)` returns an [n_year, 12] array of
area-weighted Siberian land-box means -- enough to build the JJA mean
(noise_floor), the four seasonal means (seasonal_by_run) and the month-by-month
response (monthly_lever_check) without touching the files again.

Cache invalidation follows the same rule as eval_round10_A: the signature is the
newest input mtime plus the window plus a version. BUMP CACHE_VERSION IF THE BOX,
MASK OR AVERAGING CHANGES -- a stale cache that silently disagrees with the main
table is exactly the kind of thing that has cost this campaign a round.
"""
import numpy as np, xarray as xr, os, json
from concurrent.futures import ProcessPoolExecutor

from runs import RT, LSMF, Y0, Y1

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '.box_cache')
CACHE_VERSION = 1
USE_CACHE = os.environ.get('NOCACHE', '') == ''

BOX = ((55, 75), (60, 180))
YEARS = list(range(Y0, Y1 + 1))
ACC = 3600.0

_lsm = None


def _mask():
    global _lsm
    if _lsm is None:
        _lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values
    return _lsm


def box_mean(a2d, lat, lon):
    lsm = _mask()
    ys = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
    l180 = ((lon + 180) % 360) - 180
    xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub = a2d[ii]
    m = np.isfinite(sub) & (lsm[ii] > 0.5)
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape)
    return float(np.average(sub[m], weights=w[m])) if m.any() else float('nan')


def _sig(run, var):
    newest = 0.0
    for y in YEARS:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            return None                       # incomplete -> caller skips the run
        m = os.path.getmtime(f)
        if m > newest:
            newest = m
    return f'v{CACHE_VERSION}|{var}|{Y0}-{Y1}|{newest:.0f}'


def _cached(run, var):
    """The cached array if present and current, else None. No computation."""
    sig = _sig(run, var)
    if sig is None:
        return None
    cf = os.path.join(CACHE_DIR, f'{run}__{var}.json')
    if not os.path.exists(cf):
        return None
    try:
        blob = json.load(open(cf))
        return np.array(blob['v']) if blob.get('sig') == sig else None
    except Exception:
        return None


def monthly(run, var='2t'):
    """[n_year, 12] Siberian land-box means. None if the run is incomplete.

    `2t` is returned in degC; accumulated flux variables are divided by ACC.
    """
    sig = _sig(run, var)
    if sig is None:
        return None
    os.makedirs(CACHE_DIR, exist_ok=True)
    cf = os.path.join(CACHE_DIR, f'{run}__{var}.json')
    if USE_CACHE and os.path.exists(cf):
        try:
            blob = json.load(open(cf))
            if blob.get('sig') == sig:
                return np.array(blob['v'])
        except Exception:
            pass                              # unreadable cache is not fatal
    out = []
    for y in YEARS:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        ds = xr.open_dataset(f)
        a = ds[var].values
        lat, lon = ds[var].lat.values, ds[var].lon.values
        ds.close()
        if var not in ('2t', 'tcc', 'fal', 'asn', 'rsn', 'sd'):
            a = a / ACC
        out.append([box_mean(a[m], lat, lon) for m in range(12)])
    v = np.array(out)
    if var == '2t':
        v = v - 273.15
    try:
        json.dump({'sig': sig, 'v': v.tolist()}, open(cf, 'w'))
    except Exception:
        pass
    return v


def _worker(args):
    run, var = args
    try:
        return run, monthly(run, var)
    except Exception as e:                    # one bad run must not kill the table
        print(f'  !! {run}: {type(e).__name__}: {e}')
        return run, None


def load_all(runs, var='2t', quiet=False, nproc=None):
    """[(labels, [run, n_year, 12])] for the runs that have complete output.

    Cache hits are cheap, so only the MISSES are farmed out to a process pool --
    the same pattern as eval_round10_A.py. With everything cached this returns in
    well under a second; a cold cache costs one parallel pass.
    """
    todo = [(lab, run) for lab, run in runs
            if not (USE_CACHE and _cached(run, var) is not None)]
    if todo:
        n = nproc or max(1, min(8, (os.cpu_count() or 2) - 1))
        if not quiet:
            print(f'  boxcache[{var}]: {len(todo)} to compute, '
                  f'{len(runs)-len(todo)} cached [{n} proc]')
        with ProcessPoolExecutor(max_workers=n) as ex:
            list(ex.map(_worker, [(run, var) for _, run in todo]))

    labs, data = [], []
    for lab, run in runs:
        v = monthly(run, var)
        if v is None:
            if not quiet:
                print(f'  !! {lab}: incomplete, skipped')
            continue
        labs.append(lab)
        data.append(v)
    return labs, np.array(data)
