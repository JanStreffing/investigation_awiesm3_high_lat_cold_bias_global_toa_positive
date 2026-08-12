"""Collect every stage of the ocp-tool land-sea pipeline into one npz.

Separated from the plotting so that the expensive part -- the FESOM polygon
mask, which triangulates the CORE2 mesh -- is paid once and the figure can be
re-drawn freely.

THE STAGES.  Two of the five are measured from files on disk, two are
reconstructed from code that no longer runs, and one is a proposal.  Which is
which is recorded in the npz so the figure can label them honestly:

  0  ecmwf     pristine ICMGGab45INIT, as shipped                    MEASURED
  1  legacy    the 2020-2022 rule: any cl >= 0.5 becomes land,       RECONSTRUCTED
               soil type 6.  Removed in 7dc554d.
  2  v1        current master, FESOM-driven flip + NN fill           MEASURED
  3  v2        v1 with the corrected soil map (_v2)                  MEASURED
  4  lake      v2 plus restore_flipped_lakes                         PROPOSED

THE THREE MASKS ARE NOT THE SAME MASK, which is the point of showing them
apart.  ``lsm_binary_land`` is snapshotted before the flip and ``lsm_binary_atm``
after it (ocp-tool lsm.py:352 and :422), a separation introduced in d5dcecd so
that lakes count as land for the atmosphere but not as solid land.  The FESOM
mask is a third thing again: it is what the flip is driven BY, not a product
of it.

A TRAP THIS AVOIDS.  ``gribfield_mod = gribfield[:]`` in modify_lsm is a SHALLOW
copy, so anything that reads the "pristine" list after the NN fill has run gets
the modified arrays.  Every stage here is read from its own file on disk, or
deep-copied before being altered, so no stage can contaminate another.
"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np

OCP = '/work/ab0246/a270092/software/ocp-tool'
sys.path.insert(0, OCP)
import eccodes  # noqa: E402

# Which ocean mesh. CORE3 is what the campaign actually runs; CORE2 is kept
# selectable because it is where the dateline defect lives.
MESH = os.environ.get('OCP_MESH', 'CORE3')

OUT = ('/work/ab0246/a270092/postprocessing/'
       'investigation_awiesm3_high_lat_cold_bias_global_toa_positive/'
       f'scripts/analysis/lake_pipeline_{MESH}.npz')

PRISTINE = f'{OCP}/input/openifs_input_default/ICMGGab45INIT'
V1 = f'{OCP}/output/TCO95_{MESH}/openifs_input_modified/ICMGGab45INIT_{MESH}'
V2 = f'{OCP}/output/TCO95_{MESH}/openifs_input_modified/ICMGGab45INIT_{MESH}_v2'
MASKS = f'{OCP}/output/TCO95_{MESH}/oasis_mct3_input/masks.nc'
FIELDS = ('lsm', 'slt', 'cl', 'dl')

# FLake integrates MIN(50, MAX(2, dl)); flakeene_mod.F90:211-212,238.
FLAKE_MIN, FLAKE_MAX = 2.0, 50.0
LAKE_THRESHOLD = 0.5


def read_grib(path, names):
    """shortName -> values, plus the coordinates off the lsm message."""
    out, coords = {}, None
    with open(path, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                sn = eccodes.codes_get(gid, 'shortName')
                if sn in names and sn not in out:
                    out[sn] = eccodes.codes_get_values(gid).astype('float64')
                    if sn == 'lsm' and coords is None:
                        coords = (eccodes.codes_get_array(gid, 'latitudes'),
                                  eccodes.codes_get_array(gid, 'longitudes'))
            finally:
                eccodes.codes_release(gid)
    return out, coords


print(__doc__)
print('=' * 100)
print(f'\nMESH = {MESH}')

# ------------------------------------------------------------------ measured
print('\nReading the three ICMGG files')
ecmwf, coords = read_grib(PRISTINE, FIELDS)
v1, _ = read_grib(V1, FIELDS)
v2, _ = read_grib(V2, FIELDS)
lat = np.asarray(coords[0], float)
lon = np.asarray(coords[1], float) % 360.0
n = len(ecmwf['lsm'])
for label, d in (('ecmwf', ecmwf), ('v1', v1), ('v2', v2)):
    print(f'  {label:6s} ' + '  '.join(
        f'{k}: {"ok" if k in d and len(d[k]) == n else "MISSING"}' for k in FIELDS))

# ------------------------------------------------------- reconstructed: legacy
# The rule as it stood in d5dcecd..7dc554d.  Reconstructed, not measured: the
# code that produced it is deleted, so this is what it WOULD have written to a
# pristine 48r1 file.  Soil type 6 is what the loop set; the comment beside it
# called that "SANDY CLAY LOAM", but 6 is organic in the HTESSEL table.
legacy = {k: ecmwf[k].copy() for k in FIELDS}
was_lake_ecmwf = ecmwf['cl'] >= LAKE_THRESHOLD
legacy['lsm'][was_lake_ecmwf] = 1.0
legacy['slt'][was_lake_ecmwf] = 6.0
print(f'\nLegacy rule would have promoted {int(was_lake_ecmwf.sum())} cells to '
      f'land at soil type 6')

# ------------------------------------------------------------ proposed: lake
# v2, with the lake cover and depth put back wherever the flip made a cell land
# that the pristine file had already flagged as lake.
flipped_to_land = (ecmwf['lsm'] < 0.5) & (v2['lsm'] >= 0.5)
restore = flipped_to_land & was_lake_ecmwf
lake = {k: v2[k].copy() for k in FIELDS}
lake['cl'][restore] = 1.0
lake['dl'][restore] = np.clip(ecmwf['dl'][restore], FLAKE_MIN, FLAKE_MAX)
print(f'Proposed restore touches {int(restore.sum())} cells '
      f'({int(flipped_to_land.sum())} flipped in total)')

# ------------------------------------------------------------- the FESOM mask
cache = OUT.replace('.npz', '_fesom.npy')
if os.path.exists(cache):
    fesom = np.load(cache)
    print(f'\nFESOM polygon mask from cache ({cache})')
else:
    print(f'\nBuilding the FESOM polygon mask ({MESH} mesh, slow)')
    from ocp_tool.config import load_config
    from ocp_tool.gaussian_grids import (generate_gaussian_grid,
                                         read_fesom_grid_polygon)
    cfg = load_config(f'{OCP}/configs/TCO95_{MESH}.yaml')
    grid = generate_gaussian_grid(cfg, 95)
    fesom = np.asarray(read_fesom_grid_polygon(cfg, grid, verbose=False), float)
    np.save(cache, fesom)
print(f'  FESOM mask: {int((fesom >= 0.5).sum())} land / '
      f'{int((fesom < 0.5).sum())} ocean of {len(fesom)}')

# ------------------------------------------------------------- the OASIS masks
print('\nReading the OASIS masks')
oasis = {}
try:
    import netCDF4
    with netCDF4.Dataset(MASKS) as ds:
        for key, var in (('A', 'A096.msk'), ('L', 'L096.msk'),
                         ('R', 'R096.msk'), ('land', 'TCO95-land.msk')):
            if var in ds.variables:
                oasis[key] = np.asarray(ds.variables[var][:]).ravel().astype(float)
                print(f'  {var:16s} {int((oasis[key] == 1).sum()):6d} set of '
                      f'{oasis[key].size}')
            else:
                print(f'  {var:16s} ABSENT')
except Exception as exc:                                    # noqa: BLE001
    print(f'  could not read {MASKS}: {exc}')

# ------------------------------------------------------------------------ save
save = {'lat': lat, 'lon': lon,
        'flipped_to_land': flipped_to_land, 'restore': restore,
        'fesom': fesom}
for label, d in (('ecmwf', ecmwf), ('legacy', legacy),
                 ('v1', v1), ('v2', v2), ('lake', lake)):
    for k in FIELDS:
        if k in d:
            save[f'{label}_{k}'] = d[k]
for k, v in oasis.items():
    save[f'oasis_{k}'] = v
np.savez_compressed(OUT, **save)
print(f'\nWrote {OUT}  ({len(save)} arrays, {os.path.getsize(OUT) / 1e6:.1f} MB)')
