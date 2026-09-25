"""Drake Passage volume transport, computed with tripyview's model-accurate transect path.

WHY tripyview AND NOT A HAND-ROLLED SECTION.  On an unstructured mesh a transect is not a
row of grid cells.  tripyview's do_analyse_transects works out which edges and triangles the
line actually crosses and builds a transport path from edge_cross_dxdy, so the flux it
integrates is the one FESOM's own discretisation carries.  Interpolating u,v to a straight
lon/lat line and integrating by hand does not conserve and will not agree.

The transect is the template's own Drake Passage line, 66W from 67S to 55S.

BOLUS.  This run has GM on (K_GM_max 2500), so there is an eddy-induced transport in
bolus_u/bolus_v besides the resolved one.  Both are reported.  Observational estimates are
of the total, so the bolus-included number is the one to compare: Donohue et al. 2016 give
173.3 +- 10.7 Sv, the older Cunningham et al. 2003 estimate is 134 +- 11.2 Sv.

The mesh is unrotated (abg = [0,0,0] for CORE3), so no vector rotation is applied to either
the data or the edge vectors.

NEEDS tripyview, which is an editable install whose egg-link points at a path that no longer
exists.  Run it as:
    PYTHONPATH=/work/ab0246/a270092/software/tripyview \
    /home/a/a270092/.conda/envs/esm-tools_auto_tripyview/bin/python3 \
    scripts/analysis/drake_passage_transport.py

Usage:  Y0=2090 Y1=2099 ARM=PICAL_ccnice python3 scripts/analysis/drake_passage_transport.py
"""
import os, sys, time
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
import tripyview as tpv

ROOT = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM  = os.environ.get('ARM', 'PICAL_ccnice')
Y0, Y1 = int(os.environ.get('Y0', 2090)), int(os.environ.get('Y1', 2099))
MESH = os.environ.get('MESH', '/work/ab0246/a270092/input/fesom2/core3/')
DIAG = os.path.join(MESH, 'fesom.mesh.diag.nc')
DATA = f'{ROOT}/{ARM}/outdata/fesom/'
DO_BOLUS = os.environ.get('BOLUS', '1') != '0'

TRANSECT = [[[-66, -66], [-67, -55], 'Drake Passage', False]]

print(f'{ARM}  {Y0}-{Y1}   mesh {MESH}')

# ---------------------------------------------------------------- mesh and edges
t0 = time.time()
mesh = tpv.load_mesh_fesom2(MESH, do_rot='None', do_info=False)

# fesom.mesh.diag.nc from 2.6 on names these edge_nodes / edge_face_links; tripyview's own
# template renames them to edges / edge_tri, and both are one-based on disk.
md = xr.open_dataset(DIAG, decode_cf=False)
edge        = md['edge_nodes'].values.astype(int) - 1
edge_tri    = md['edge_face_links'].values.astype(int) - 1
edge_dxdy   = md['edge_cross_dxdy'].values
edge_dxdy_l = np.array([edge_dxdy[0, :], edge_dxdy[1, :]])
edge_dxdy_r = np.array([edge_dxdy[2, :], edge_dxdy[3, :]])
edge_dxdy_r[:, edge_tri[1, :] < 0] = 0.0      # boundary edge has no right triangle
md.close()
print(f'  mesh + edges loaded, {edge.shape[1]} edges, {time.time()-t0:.0f}s')

transects = tpv.do_analyse_transects(TRANSECT, mesh, edge, edge_tri,
                                     edge_dxdy_l, edge_dxdy_r, do_rot=False, do_info=False)
print(f'  transect path: {len(transects)} transect(s), '
      f'{len(transects[0]["edge_cut_i"])} crossed edges')

# ---------------------------------------------------------------- transport
def load_uv(nu, nv):
    """Build the (time, elem, nz1) dataset calc_transect_Xtransp expects.

    tripyview's own load_data_fesom2 needs xr.coders, which arrived after the xarray in
    this environment, so the two files are opened directly instead. The output files call
    the vertical dimension nz with 47 entries; that is FESOM's mid-layer count, which the
    mesh diagnostics call nz1, and tripyview keys its layer thickness off the name.
    """
    fu = [f'{DATA}/{nu}.fesom.{y}.nc' for y in range(Y0, Y1 + 1)]
    fv = [f'{DATA}/{nv}.fesom.{y}.nc' for y in range(Y0, Y1 + 1)]
    if not all(os.path.exists(f) for f in fu + fv):
        return None
    ku = dict(combine='by_coords', parallel=False, chunks={'time': 1})
    u = xr.open_mfdataset(fu, **ku)[nu]
    v = xr.open_mfdataset(fv, **ku)[nv]
    d = xr.Dataset({nu: u, nv: v})
    if 'nz' in d.dims and d.sizes['nz'] == 47:
        d = d.rename({'nz': 'nz1'})
    # calc_transect_Xtransp reads data_uv['lon'] / ['lat'] for the element positions
    d = d.drop_vars([c for c in d.coords if c not in ('time', 'lon', 'lat')], errors='ignore')
    # tripyview indexes the velocities as [depth, along-transect] and scales by dz along
    # axis 0, so the vertical dimension has to come before the horizontal one.
    d = d.transpose('time', 'nz1', 'elem')
    d.attrs['descript'] = ARM
    return d

def transport(nu, nv):
    d = load_uv(nu, nv)
    if d is None:
        return None
    cs = tpv.calc_transect_Xtransp(mesh, d, transects, do_rot=False,
                                   do_tarithm='mean', do_info=False, client=None)
    return cs[0]

res = {}
for lab, nu, nv in (('resolved', 'u', 'v'),) + ((('bolus', 'bolus_u', 'bolus_v'),) if DO_BOLUS else ()):
    t = time.time()
    cs = transport(nu, nv)
    if cs is None:
        print(f'  {lab}: {nu}/{nv} not available'); continue
    key = list(cs.keys())[0]
    a = cs[key]
    res[lab] = a
    print(f'  {lab:8s} computed from {nu}/{nv} in {time.time()-t:.0f}s  '
          f'(var {key}, dims {dict(a.sizes)})')

if not res:
    sys.exit('no transport computed')

# ---------------------------------------------------------------- report
# calc_transect_Xtransp returns Vflx already in Sv, one value per (level, path segment).
# Do not rescale it.
assert all(a.attrs.get('units') == 'Sv' for a in res.values()), 'unexpected units'
SV = 1.0
print()
print(f'{"component":<12}{"net":>10}{"eastward":>12}{"westward":>12}   [Sv]')
tot = None
for lab, a in res.items():
    v = np.asarray(a.values, float)
    net  = np.nansum(v) * SV
    east = np.nansum(np.where(v > 0, v, 0.0)) * SV
    west = np.nansum(np.where(v < 0, v, 0.0)) * SV
    tot = v if tot is None else tot + v
    print(f'{lab:<12}{net:>10.1f}{east:>12.1f}{west:>12.1f}')
if len(res) > 1:
    print(f'{"TOTAL":<12}{np.nansum(tot)*SV:>10.1f}'
          f'{np.nansum(np.where(tot>0,tot,0.0))*SV:>12.1f}'
          f'{np.nansum(np.where(tot<0,tot,0.0))*SV:>12.1f}')
print()
print('  observed: Donohue 2016  173.3 +- 10.7 Sv;  Cunningham 2003  134 +- 11.2 Sv')
