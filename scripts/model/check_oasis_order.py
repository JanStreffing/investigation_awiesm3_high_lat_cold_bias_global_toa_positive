"""Check an OASIS set's FESOM order against a FESOM partition (the sufficient test).

OASIS 'feom' order = concatenation over ranks r of the first myDim entries of
dist_N/my_listRRRRR.out (line 2 = myDim, the list starts at token 4).  grids.nc
feom.lon/lat must then equal nod2d.out coordinates in that order; 1609 polar nodes
differ by < 0.006 deg from rotation roundoff, a wrong order is off by up to 180 deg.
Also runs the README's sorted-hash test (nothing lost) on the restart files.

Usage: check_oasis_order.py <mesh_dir> <nproc> <oasis_out_dir> <oasis_in_dir>
"""
import sys, hashlib, numpy as np
from netCDF4 import Dataset
mesh, npr, out, inp = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
with open(f'{mesh}/nod2d.out') as f:
    n = int(f.readline())
    xy = np.loadtxt(f, max_rows=n, usecols=(1, 2))
order = []
for r in range(npr):
    tok = open(f'{mesh}/dist_{npr}/my_list{r:05d}.out').read().split()
    my = int(tok[1])
    order.extend(int(t) - 1 for t in tok[3:3 + my])
order = np.array(order)
assert len(order) == n and len(np.unique(order)) == n, (len(order), n)
with Dataset(f'{out}/grids.nc') as g:
    lon = np.asarray(g['feom.lon'][:]).ravel(); lat = np.asarray(g['feom.lat'][:]).ravel()
dlon = np.abs(((lon - xy[order, 0]) + 180) % 360 - 180); dlat = np.abs(lat - xy[order, 1])
bad = (dlon > 0.01) | (dlat > 0.01)
print(f'nodes {n}; max |dlon| {dlon.max():.4f}  max |dlat| {dlat.max():.4f}  off by > 0.01 deg: {bad.sum()}')
def vh(ds, vn):
    a = np.asarray(ds.variables[vn][...]).ravel()
    return hashlib.sha256(np.sort(a).tobytes()).hexdigest()[:16]
for fn in ('rstas.nc', 'rstos.nc', 'vegin.nc'):
    a, b = Dataset(f'{inp}/{fn}'), Dataset(f'{out}/{fn}')
    vs = sorted(set(a.variables) & set(b.variables))
    ok = sum(vh(a, v) == vh(b, v) for v in vs)
    print(f'{fn}: sorted-hash equal {ok}/{len(vs)}')
