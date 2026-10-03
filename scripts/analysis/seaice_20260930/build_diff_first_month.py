"""Is the new FESOM build numerically the same as the old one? First-month field differences.

The 766020ec -> 53e095ec diff is meant to be refactors (lock removal, node-owned gathers,
EVP halo overlap, IDEMIX2 gathers).  An exact refactor gives first-month differences no
larger than chaotic divergence; a real change gives larger, systematic ones, in particular in
the deep diffusivity that IDEMIX sets.  Reference for chaos: two runs of the SAME build and
time step that differ only in OpenMP thread count (nolockB1 vs nolockB4).
Metric per depth band: area-weighted mean |a-b| / mean(|a|+|b|)/2 of the January mean, and
the mean signed log10 ratio (systematic bias), for Kv, temp, salt.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4'
PAIRS = [('noise: same new build 1800 s, 1 vs 4 threads', ('TEST_nolockB1', 2100), ('TEST_nolockB4', 2100)),
         ('STEP, old build: perfE 1200 s vs perfF 1800 s', ('TEST_perfE', 2100), ('TEST_perfF', 2100)),
         ('STEP, old build: ob1200 vs ob1800', ('PICAL_crunveg_ob1200', 2120), ('PICAL_crunveg_ob1800', 2120)),
         ('STEP, new build: ts1200 vs ihf0', ('PICAL_crunveg_ts1200', 2120), ('PICAL_crunveg_ihf0', 2120)),
         ('build at 1800 s: old perfF vs new nolockB1', ('TEST_perfF', 2100), ('TEST_nolockB1', 2100)),
         ('build at 1200 s: old ob1200 vs new ts1200', ('PICAL_crunveg_ob1200', 2120), ('PICAL_crunveg_ts1200', 2120))]


def jan(exp, y, v):
    with xr.open_dataset(f'{R}/{exp}/outdata/fesom/{v}.fesom.{y}.nc', decode_times=False) as d:
        a = d[v].isel(time=0).values.astype('f8')
        dep = d[[c for c in d[v].dims if c.startswith('nz')][0]].values
    return np.where(np.abs(a) > 1e30, np.nan, a), dep


for lab, (ea, ya), (eb, yb) in PAIRS:
    print(lab)
    for v in ('Kv', 'Av', 'temp', 'salt'):
        try:
            a, dep = jan(ea, ya, v); b, _ = jan(eb, yb, v)
        except Exception as e:
            print(f'   {v}: {e}'); continue
        if a.shape[0] == len(dep): a, b = a.T, b.T          # -> (node, depth)
        out = []
        for z0, z1 in ((0, 100), (100, 500), (500, 1500), (1500, 7000)):
            kz = (np.abs(dep) >= z0) & (np.abs(dep) < z1)
            x, y = a[:, kz], b[:, kz]
            ok = np.isfinite(x) & np.isfinite(y)
            if v in ('Kv', 'Av'):
                ok &= (x > 0) & (y > 0)
                rel = np.mean(np.abs(x[ok] - y[ok])) / np.mean(0.5 * (x[ok] + y[ok]))
                bias = np.mean(np.log10(y[ok] / x[ok]))
                out.append(f'{z0}-{z1}m rel {rel:6.3f} log10(b/a) {bias:+7.4f}')
            else:
                out.append(f'{z0}-{z1}m mean|d| {np.mean(np.abs(x[ok]-y[ok])):8.5f} mean d {np.mean(y[ok]-x[ok]):+9.5f}')
        print(f'   {v:5s} ' + '   '.join(out), flush=True)
