"""Is there a Southern Ocean SUBSURFACE warm bias, and could it be melting the ice?

WHY THIS IS OPEN AT ALL.  release_evaluation_tool2 produces 0m/100m/1000m/4000m
temperature-vs-PHC3 figures, and for every run in this campaign they are empty -- grey
maps reading rmsd=nan, bias=nan, byte-identical across 11E/11G/11N/11P/11R.  So the
ocean interior has never actually been validated here, and the question of whether the
SH sea ice is melted from above (cloud shortwave) or from below (ocean heat) has never
been decidable.

WHY IT MATTERS NOW.  The melt-season energy excess at 90-60S is +6.98 W/m2 of cloud
shortwave, which argues melt-from-above.  But the MAM warm anomaly is +4.48 W/m2 of
LONGWAVE with almost no shortwave, in the refreeze season when the sun is leaving --
ice that will not refreeze in autumn is the classic signature of heat arriving from
underneath.  Those want different levers, so the interior has to be measured.

METHOD, and why no interpolation is needed.  The PHC3 climatology staged under obs/phc3
is itself in FESOM format, but on a DIFFERENT mesh: 126858 nodes against the CORE3 mesh's
211567.  The two cannot be differenced node-by-node.  They do however share an identical
47-level vertical grid (2.5, 7.5, 15, 25 ... m), so both are binned onto a common
(latitude, depth) grid and compared as zonal means.  That sidesteps mesh interpolation
entirely, at the cost of saying nothing about horizontal structure within a band.

Node areas come from each mesh's own fesom.mesh.diag.nc (nod_area), so the zonal means
are area-weighted rather than node-count-weighted -- the CORE3 mesh refines towards the
coast, and an unweighted mean would be dominated by the shelf.

PHC3 node coordinates live in the mesh diag as `nodes` (2, nod2) in RADIANS, not in the
data file.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
PHC = '/work/ab0246/a270092/obs/phc3'
MESH = '/work/ab0246/a270092/input/fesom2/core3/fesom.mesh.diag.nc'
Y0, Y1 = 1380, 1389
ARMS = [('11P', f'{R}/11P'), ('11R', f'{R}/11R'),
        ('11N', f'{R}/11N'), ('11Q', f'{R}/11Q')]
BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45)]
DEPTHS = [0, 50, 100, 200, 500, 1000, 2000]


def node_meta(diag):
    """The two meshes store node latitude differently: the PHC3 diag has `nodes`
    (2, nod2) in RADIANS, the CORE3 diag has explicit `lat`/`lon` in degrees.  Handle
    both rather than assume, and sanity-check the range afterwards."""
    with xr.open_dataset(diag, decode_times=False) as d:
        if 'lat' in d:
            lat = np.asarray(d['lat'].values)
        else:
            lat = np.asarray(d['nodes'].values)[1] * (180.0 / np.pi)
        if np.nanmax(np.abs(lat)) < 3.2:             # still radians
            lat = lat * (180.0 / np.pi)
        area = np.asarray(d['nod_area'].values)      # (nz, nod2)
    if not (-91 < np.nanmin(lat) and np.nanmax(lat) < 91):
        raise SystemExit(f'{diag}: latitude range {np.nanmin(lat)}..{np.nanmax(lat)}')
    return lat, area


def zonal(temp, lat, area, nz_vals):
    """temp (nod2, nz) -> dict[(band, depth)] = area-weighted mean."""
    out = {}
    for bname, lo, hi in BANDS:
        sel = (lat >= lo) & (lat < hi)
        for dtgt in DEPTHS:
            iz = int(np.argmin(np.abs(nz_vals - dtgt)))
            w = area[min(iz, area.shape[0] - 1)][sel]
            v = temp[sel, iz]
            k = np.isfinite(v) & np.isfinite(w) & (w > 0) & (v > -5) & (v < 40)
            out[(bname, dtgt)] = float(np.average(v[k], weights=w[k])) if k.any() else np.nan
    return out


# --- observed ---
with xr.open_dataset(f'{PHC}/temp.fesom.1958.nc', decode_times=False) as d:
    o = np.squeeze(d['temp'].values)                 # (nod2, nz1)
    nzv = np.asarray(d['nz1'].values)
olat, oarea = node_meta(f'{PHC}/fesom.mesh.diag.nc')
obs = zonal(o, olat, oarea, nzv)

mlat, marea = node_meta(MESH)

print(__doc__)
print('=' * 92)
print(f'\nSouthern Ocean temperature vs PHC3, {Y0}-{Y1}, area-weighted zonal means [degC]')
for bname, _, _ in BANDS:
    print(f'\n--- {bname} ---')
    print(f'{"depth":>7} {"PHC3":>8} | ' + ' | '.join(f'{a:>14}' for a, _ in ARMS))
    rows = {}
    for lab, path in ARMS:
        acc = []
        for y in range(Y0, Y1 + 1):
            p = f'{path}/outdata/fesom/temp.fesom.{y}.nc'
            if not os.path.exists(p):
                acc = None
                break
            with xr.open_dataset(p, decode_times=False) as d:
                acc.append(np.asarray(d['temp'].values).mean(0))     # annual, (nod2, nz)
        rows[lab] = zonal(np.mean(acc, axis=0), mlat, marea, nzv) if acc else None
    for dtgt in DEPTHS:
        line = f'{dtgt:>7} {obs[(bname, dtgt)]:8.2f} | '
        cells = []
        for lab, _ in ARMS:
            r = rows[lab]
            cells.append('           n/a' if r is None
                         else f'{r[(bname, dtgt)]:7.2f} ({r[(bname, dtgt)] - obs[(bname, dtgt)]:+5.2f})')
        print(line + ' | '.join(cells))
