"""Where does the heat GM 2500 keeps out of mid-depth enter, where is it stored, and by which
route does it get there? Annual heat budget by basin and depth for one run (CORE3 mesh).

Regions (pyfesom2 MOCBasins.geojson, plus latitude):
  SO      south of 34S
  SATL    Atlantic_MOC, 34S-0
  NATL    Atlantic_MOC, 0-45N
  SPNA    Atlantic_MOC, 45-66N
  NORD    north of 66N (Nordic Seas + Arctic, any basin)
  IP      IndoPacific_MOC north of 34S (and south of 66N)
  OTHER   everything else (Mediterranean, marginal seas outside the polygons)
Output rows: exp,year,kind,region,band,value
  ohc   [J]   heat content 4.1e6*T*V per region and depth band (0-700, 700-2000, 2000+ m)
  sfc   [W]   heat entering through the surface, -fh (fh is positive upward)
  tend / adv / diff [W]  column-integrated CMIP budget terms (opottemptend, opottemprmadvect,
        opottempdiff, W m-2 per column) summed over the region. With the linear free surface the
        column advection also carries a surface term (its global sum is not 0), so use
        tend - sfc as the lateral import into a region (global sum about 0).
  mht_res / mht_bol [W]  northward heat transport across a latitude in the Atlantic (Atlantic_MOC)
        per depth band: elements whose centroid lies within
        +-1 deg of the latitude, sum(4.1e6 * vT * A_elem * dz) / band width. Resolved from vtemp;
        bolus from bolus_v * T (annual means multiplied, so the eddy part of the bolus term is
        left out). region field = '<basin>@<lat>'.
Usage: python heat_pathways.py <mesh.diag.nc> <MOCBasins.geojson> <outdata_fesom> <exp> <y0> <y1> <out.csv>
y0 is the year whose heat content anchors the first change; budget terms are written from y0+1.
"""
import sys, json, numpy as np
from netCDF4 import Dataset
from matplotlib.path import Path
diag, geoj, d, exp, y0, y1, outf = sys.argv[1:8]; y0, y1 = int(y0), int(y1)
RHOCP = 4.1e6
with Dataset(diag) as nc:
    narea = np.asarray(nc.variables['nod_area'][:], 'f8')            # (48, nod2)
    earea = np.asarray(nc.variables['elem_area'][:], 'f8')           # (elem)
    fnod = np.asarray(nc.variables['face_nodes'][:], 'i8')           # (3, elem)
    zb = np.abs(np.asarray(nc.variables['nz'][:], 'f8'))
    nlon = np.asarray(nc.variables['lon'][:], 'f8'); nlat = np.asarray(nc.variables['lat'][:], 'f8')
if np.abs(nlon).max() < 7: nlon, nlat = np.rad2deg(nlon), np.rad2deg(nlat)
nlon = ((nlon + 180) % 360) - 180
if fnod.min() == 1: fnod -= 1
with Dataset(f'{d}/vtemp.fesom.{y0+1}.nc') as nc:
    elon = ((np.asarray(nc.variables['lon'][:], 'f8') + 180) % 360) - 180; elat = np.asarray(nc.variables['lat'][:], 'f8')
zmid = 0.5 * (zb[:-1] + zb[1:])
BAND = {'0-700': zmid < 700, '700-2000': (zmid >= 700) & (zmid < 2000), '2000+': zmid >= 2000, 'all': zmid >= 0}
G = {f['properties']['name']: f['geometry'] for f in json.load(open(geoj))['features']}
def inpoly(name, lon, lat):
    # exterior rings only (the basin polygons have no holes that matter here)
    g = G[name]; polys = [g['coordinates']] if g['type'] == 'Polygon' else g['coordinates']
    xy = np.column_stack([lon, lat]); m = np.zeros(len(lon), bool)
    for p in polys: m |= Path(np.asarray(p[0])).contains_points(xy)
    return m
natl, nip = inpoly('Atlantic_MOC', nlon, nlat), inpoly('IndoPacific_MOC', nlon, nlat)
eatl, eip = inpoly('Atlantic_MOC', elon, elat), inpoly('IndoPacific_MOC', elon, elat)
REG = {'SO': nlat < -34, 'NORD': nlat >= 66,
       'SATL': natl & (nlat >= -34) & (nlat < 0), 'NATL': natl & (nlat >= 0) & (nlat < 45),
       'SPNA': natl & (nlat >= 45) & (nlat < 66), 'IP': nip & (nlat >= -34) & (nlat < 66)}
used = np.zeros_like(nlat, bool)
for k in ('SO', 'NORD', 'SATL', 'NATL', 'SPNA', 'IP'):
    REG[k] = REG[k] & ~used; used |= REG[k]
REG['OTHER'] = ~used
# Only the Atlantic is kept: Indo-Pacific and global sections carry net throughflow (ITF, Bering),
# so v*T against 0 degC is not a heat transport there, and the equatorial band is meaningless.
EB = {'ATL': eatl}
LATS = (-34, -20, 20, 30, 45, 55, 66)
# layer thickness: time mean hnode of the anchor year at nodes, element = mean of its 3 nodes
with Dataset(f'{d}/hnode.fesom.{y0}.nc') as nc:
    h = np.asarray(nc.variables['hnode'][:], 'f8').mean(0)
h = h if h.shape[0] == 47 else h.T
h = np.where(np.isfinite(h) & (np.abs(h) < 1e5), h, 0.0)
vol = h * narea[:47]
he = h[:, fnod].mean(1)                                               # (47, elem)
def annual(v, f):
    with Dataset(f) as nc:
        x = nc.variables[v]; acc = None
        for m in range(x.shape[0]):
            a = np.asarray(x[m], 'f8'); a[~np.isfinite(a) | (np.abs(a) > 1e10)] = 0.0
            acc = a if acc is None else acc + a
        return acc / x.shape[0]
out = open(outf, 'a')
def w(year, kind, reg, band, val): out.write(f'{exp},{year},{kind},{reg},{band},{val:.6e}\n')
for y in range(y0, y1 + 1):
    t = annual('temp', f'{d}/temp.fesom.{y}.nc'); t = t if t.shape[0] == 47 else t.T
    hc = RHOCP * t * vol
    for r, m in REG.items():
        for b, k in BAND.items(): w(y, 'ohc', r, b, hc[k][:, m].sum())
    if y == y0: out.flush(); continue
    for v, kind in (('fh', 'sfc'), ('opottemptend', 'tend'), ('opottemprmadvect', 'adv'), ('opottempdiff', 'diff')):
        a = annual(v, f'{d}/{v}.fesom.{y}.nc')
        if kind == 'sfc': a = -a
        for r, m in REG.items(): w(y, kind, r, 'col', (a[m] * narea[0][m]).sum())
    te = t[:, fnod].mean(1)                                            # T at elements (47, elem)
    vt = annual('vtemp', f'{d}/vtemp.fesom.{y}.nc'); vt = vt if vt.shape[0] == 47 else vt.T
    vb = annual('bolus_v', f'{d}/bolus_v.fesom.{y}.nc'); vb = vb if vb.shape[0] == 47 else vb.T
    for la in LATS:
        sel = np.abs(elat - la) < 1.0
        width = 2.0 * 111.2e3
        for bn, bm in EB.items():
            s = sel & bm
            if not s.any(): continue
            res = RHOCP * vt[:, s] * he[:, s] * earea[s] / width
            bol = RHOCP * vb[:, s] * te[:, s] * he[:, s] * earea[s] / width
            for b, k in BAND.items():
                w(y, 'mht_res', f'{bn}@{la}', b, res[k].sum()); w(y, 'mht_bol', f'{bn}@{la}', b, bol[k].sum())
    out.flush(); print(f'{exp} {y} ok', flush=True)
