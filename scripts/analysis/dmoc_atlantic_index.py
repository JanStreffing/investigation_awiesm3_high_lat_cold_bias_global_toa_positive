"""Atlantic overturning in DENSITY space (sigma2) from FESOM's online density-class diagnostics.

WHY.  The z-space AMOC at 26.5N fell by ~2 Sv at 2120 and the hemispheric-GM run sits 1.1 Sv below
its same-years control, although its Labrador Sea convects.  A z-space streamfunction mixes the
overturning with the gyre where isopycnals slope across the basin (north of ~40N), and it does not
see the eddy-induced (GM bolus) transport.  The density-space streamfunction does both.

METHOD.  std_dens_DIV(node, sigma2 class) is the divergence of the volume transport within each
density class, binned online every step, so the annual mean keeps the seasonal and transient
correlations.  For the basin north of latitude y (Atlantic_MOC mask of pyfesom2: Atlantic, Arctic,
Nordic Seas), the northward transport across y in class k is the sum of the divergence over all
nodes north of y.  Cumulating over classes from the densest upwards gives psi(y, sigma).
std_dens_DIVbolus is the same for the GM eddy-induced velocity; residual = Eulerian + bolus.
The sign is chosen so that the upper (northward light water, southward dense water) cell is
positive.  Only the two divergence terms are used; the surface transformation and the dV/dt terms
of the full framework are not needed for the transport itself.

Per run and year, cached: maximum of psi over density at 26.5N, 45N and 55N, and the density of the
maximum, for the Eulerian and the residual flow.

Usage: dmoc_atlantic_index.py <tag> <exp root> <y0> <y1>     (reval environment, one process per call)
       dmoc_atlantic_index.py --table                        (print decadal means of everything cached)
"""
import os
import sys
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
CACHE = f'{P}/reval/spinup_cache/dmoc'
MESH = f'{P}/reval_obs/mesh/core3/'
LATS = np.arange(-30.0, 80.5, 1.0)
PICK = (26.5, 45.0, 55.0)
SIG_MIN = 35.0          # ignore the light surface classes when looking for the overturning maximum


def psi_from_div(div, lat, mask):
    """div (node, class) [m3/s] -> psi (lat, class) [Sv], cumulated from the densest class up."""
    out = np.zeros((len(LATS), div.shape[1]))
    for i, y in enumerate(LATS):
        m = mask & (lat >= y)
        out[i] = div[m].sum(axis=0)
    return np.cumsum(out[:, ::-1], axis=1)[:, ::-1] * 1e-6


def year_index(root, yr, lat, mask):
    import xarray as xr
    fields = {}
    for name in ('std_dens_DIV', 'std_dens_DIVbolus'):
        fn = f'{root}/outdata/fesom/{name}.fesom.{yr}.nc'
        if not os.path.exists(fn):
            return None
        with xr.open_dataset(fn) as ds:
            sig = ds['std_dens'].values
            fields[name] = np.nan_to_num(ds[name].mean('time').transpose('nod2', 'std_dens').values.astype(np.float64))
    eul = psi_from_div(fields['std_dens_DIV'], lat, mask)
    tot = psi_from_div(fields['std_dens_DIV'] + fields['std_dens_DIVbolus'], lat, mask)
    return sig, eul, tot


if len(sys.argv) > 1 and sys.argv[1] != '--table':
    import pyfesom2 as pf
    from pyfesom2.ut import get_mask
    tag, root, y0, y1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    os.makedirs(f'{CACHE}/{tag}', exist_ok=True)
    mesh = pf.load_mesh(MESH)
    mask = get_mask(mesh, 'Atlantic_MOC')
    for yr in range(y0, y1 + 1):
        fn = f'{CACHE}/{tag}/dmoc_{yr}.npz'
        if os.path.exists(fn):
            continue
        r = year_index(root, yr, mesh.y2, mask)
        if r is None:
            print(f'{tag} {yr}: no density-class output', flush=True)
            continue
        sig, eul, tot = r
        np.savez_compressed(fn, sig=sig, lat=LATS, eul=eul.astype(np.float32), res=tot.astype(np.float32))
        j = int(np.argmin(np.abs(LATS - 26.5)))
        k = sig >= SIG_MIN
        sgn = 1.0 if abs(eul[j, k].max()) >= abs(eul[j, k].min()) else -1.0
        print(f'{tag} {yr}: 26.5N Eulerian {sgn * (sgn * eul[j, k]).max():+.2f} Sv, residual {sgn * (sgn * tot[j, k]).max():+.2f} Sv', flush=True)
    sys.exit(0)

# ------------------------------------------------------------------ table
import glob
import pandas as pd
rows = []
for d in sorted(glob.glob(f'{CACHE}/*')):
    tag = os.path.basename(d)
    for fn in sorted(glob.glob(f'{d}/dmoc_*.npz')):
        z = np.load(fn)
        sig, k = z['sig'], z['sig'] >= SIG_MIN
        row = dict(run=tag, year=int(fn[-8:-4]))
        for nm in ('eul', 'res'):
            psi = z[nm]
            j26 = int(np.argmin(np.abs(LATS - 26.5)))
            sgn = 1.0 if abs(psi[j26, k].max()) >= abs(psi[j26, k].min()) else -1.0
            for y in PICK:
                j = int(np.argmin(np.abs(LATS - y)))
                col = sgn * psi[j, k]
                row[f'{nm}_{y:g}N'] = col.max()
                row[f'{nm}_sig_{y:g}N'] = sig[k][int(col.argmax())]
            band = (LATS >= 30) & (LATS <= 60)
            row[f'{nm}_max30-60N'] = (sgn * psi[band][:, k]).max()
        rows.append(row)
df = pd.DataFrame(rows)
df.to_csv(f'{REPO}/data/clim/dmoc_atlantic_index_annual.csv', index=False, float_format='%.4g')
cols = ['eul_26.5N', 'res_26.5N', 'eul_45N', 'res_45N', 'eul_55N', 'res_55N', 'res_max30-60N', 'res_sig_45N']
out = ['Atlantic overturning in sigma2 space, maximum over density classes >= 35.0 [Sv]; eul = Eulerian, res = Eulerian + GM bolus.',
       'Decadal means; the last column is the density of the 45N residual maximum.',
       f"{'run':22s} {'years':10s} {'n':>3s} " + ' '.join(f'{c:>13s}' for c in cols)]
for tag, g in df.groupby('run', sort=False):
    y0 = int(g.year.min())
    for a in range(y0, int(g.year.max()) + 1, 10):
        s = g[(g.year >= a) & (g.year < a + 10)]
        out.append(f'{tag:22s} {a}-{int(s.year.max())}  {len(s):3d} ' + ' '.join(f'{s[c].mean():13.2f}' for c in cols))
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/dmoc_atlantic_index.txt', 'w').write(txt + '\n')
