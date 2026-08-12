"""What is actually at the cells ocp-tool flips from ocean to land?

WHY THIS EXISTS.  ocp-tool reconciles the OpenIFS land-sea mask with the FESOM
mesh.  Where IFS says ocean and FESOM has no wet node, the cell is made land and
its whole surface column is rebuilt from the nearest stable LAND neighbour
(``fill_flipped_from_nearest_neighbour``).  Lake cover ``cl`` is read by
``read_grib_fields`` and never written, so it is carried along by the NN fill
like any other field -- it inherits the neighbour's value, which for a dry
inland neighbour is zero.

The question is whether those cells would be better represented as FLake lakes
than as dry soil.  Answering it needs three measurements that no amount of
reading the source can supply:

  1. What ``cl`` and ``dl`` hold in the PRISTINE input at the flipped cells.
     If the input already calls them lake, the flip is destroying information
     that is right there in the file.
  2. What they hold AFTER the flip, i.e. what the NN fill did to them.
  3. Whether ``dl`` is usable: FLake clips the depth it integrates to
     [2, 50] m (``flakeene_mod.F90:211-212, 238``), so any value outside that
     band is equivalent to the clip, and a zero or missing value silently
     becomes 2 m rather than failing.

A TRAP THIS AVOIDS.  The reduced Gaussian grid is not a lat-lon array; field
index i means nothing without the grid description.  The coordinates are taken
from the GRIB message itself (``latitudes``/``longitudes``), which is by
construction the same ordering as the values array.
"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np
import eccodes

PRISTINE = ('/work/ab0246/a270092/software/ocp-tool/input/'
            'openifs_input_default/ICMGGab45INIT')
FLIPPED = ('/work/ab0246/a270092/software/ocp-tool/output/TCO95_CORE2/'
           'openifs_input_modified/ICMGGab45INIT_CORE2_v2')
WANT = ('lsm', 'cl', 'dl', 'slt', 'cvh', 'cvl', 'lmlt', 'lblt', 'licd', 'sst')

# FLake integrates a clipped depth; see flakeene_mod.F90.
FLAKE_D_MIN, FLAKE_D_MAX = 2.0, 50.0

REGIONS = {
    'Caspian/Aral':  (35.0, 50.0,  46.0,  62.0),
    'Great Lakes':   (41.0, 50.0, 265.0, 285.0),
    'Baltic':        (53.0, 66.0,  10.0,  30.0),
    'Hudson Bay':    (51.0, 65.0, 265.0, 290.0),
    'Arctic coast':  (66.0, 90.0,   0.0, 360.0),
    'Antarctic':    (-90.0, -60.0,  0.0, 360.0),
    'Tropics':      (-23.5,  23.5,  0.0, 360.0),
}


def read(path, names, with_coords=False):
    """shortName -> values, for a reduced-Gaussian ICMGG."""
    out, coords = {}, None
    with open(path, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                sn = eccodes.codes_get(gid, 'shortName')
                if sn in names and sn not in out:
                    out[sn] = eccodes.codes_get_values(gid)
                    if with_coords and coords is None and sn == 'lsm':
                        coords = (eccodes.codes_get_array(gid, 'latitudes'),
                                  eccodes.codes_get_array(gid, 'longitudes'))
            finally:
                eccodes.codes_release(gid)
    return (out, coords) if with_coords else out


print(__doc__)
print('=' * 100)

(a, coords), b = read(PRISTINE, WANT, with_coords=True), read(FLIPPED, WANT)
missing = [n for n in WANT if n not in a or n not in b]
if 'lsm' in missing or 'cl' in missing or 'dl' in missing:
    sys.exit(f'  required field(s) absent: {missing}')
if missing:
    print(f'\n  (absent, skipped: {missing})')

n = len(a['lsm'])
lat = np.asarray(coords[0], float)[:n]
lon = np.asarray(coords[1], float)[:n] % 360.0

flip = np.where((a['lsm'] < 0.5) & (b['lsm'] >= 0.5))[0]   # ocean -> land
back = np.where((a['lsm'] >= 0.5) & (b['lsm'] < 0.5))[0]   # land  -> ocean
print(f'\n  {n} cells; {len(flip)} flipped ocean->land, {len(back)} land->ocean\n')
if len(flip) == 0:
    sys.exit('  nothing flipped to land -- wrong file pair?')

# ------------------------------------------------------------ 1. before / after
print('=' * 100)
print('\nAT THE FLIPPED CELLS: pristine input vs what ocp-tool wrote\n')
print(f'  {"field":6s} {"in: mean":>10s} {"in: max":>10s} {"in: >0.5":>9s}'
      f'   {"out: mean":>10s} {"out: max":>10s} {"out: >0.5":>9s}')
for name in ('cl', 'dl', 'slt', 'cvh', 'cvl'):
    if name not in a:
        continue
    x, y = a[name][flip], b[name][flip]
    print(f'  {name:6s} {x.mean():10.3f} {x.max():10.2f} {(x > 0.5).sum():9d}'
          f'   {y.mean():10.3f} {y.max():10.2f} {(y > 0.5).sum():9d}')

lake_in = a['cl'][flip] >= 0.5
print(f'\n  Of {len(flip)} flipped cells, {lake_in.sum()} were ALREADY lake '
      f'(cl >= 0.5) in the pristine input.')
print(f'  After the flip, {(b["cl"][flip] >= 0.5).sum()} still are.')

# ------------------------------------------------------------ 2. is dl usable
print('\n' + '=' * 100)
print(f'\nIS dl USABLE?  FLake clips to [{FLAKE_D_MIN}, {FLAKE_D_MAX}] m, so a value '
      f'outside\nthat band is equivalent to the clip and a zero becomes {FLAKE_D_MIN} m.\n')
d_in = a['dl'][flip]
for label, sel in (('== 0 (would clip up to the minimum)', d_in <= 0.0),
                   (f'in [0, {FLAKE_D_MIN})', (d_in > 0) & (d_in < FLAKE_D_MIN)),
                   (f'in [{FLAKE_D_MIN}, {FLAKE_D_MAX}]  (used as given)',
                    (d_in >= FLAKE_D_MIN) & (d_in <= FLAKE_D_MAX)),
                   (f'> {FLAKE_D_MAX} (clipped down)', d_in > FLAKE_D_MAX)):
    print(f'  {label:48s} {int(sel.sum()):5d}')
glob = a['dl']
wet = glob[a['lsm'] < 0.5]
print(f'\n  dl over the whole globe:  min {glob.min():.2f}  max {glob.max():.2f}  '
      f'median {np.median(glob):.2f}')
print(f'  dl at pristine OCEAN cells: median {np.median(wet):.2f}  '
      f'-> the field carries bathymetry, not just lake depth')
truelake = a['dl'][a['cl'] >= 0.5]
if truelake.size:
    print(f'  dl where cl >= 0.5 (real lakes): min {truelake.min():.2f}  '
          f'median {np.median(truelake):.2f}  max {truelake.max():.2f}')
    inband = ((truelake >= FLAKE_D_MIN) & (truelake <= FLAKE_D_MAX)).mean()
    print(f'  ...of which {100 * inband:.0f} % already sit inside FLake\'s band')

# ------------------------------------------------------------ 3. geography
print('\n' + '=' * 100)
print('\nWHERE THE FLIPPED CELLS ARE, and what the input called them\n')
print(f'  {"region":16s} {"cells":>6s} {"cl>=0.5 in":>11s} {"dl in-band":>11s} '
      f'{"out slt":>9s}')
claimed = np.zeros(len(flip), bool)
for name, (la0, la1, lo0, lo1) in REGIONS.items():
    m = ((lat[flip] >= la0) & (lat[flip] < la1)
         & (lon[flip] >= lo0) & (lon[flip] < lo1) & ~claimed)
    if not m.any():
        continue
    claimed |= m
    idx = flip[m]
    d = a['dl'][idx]
    slt_out = b['slt'][idx] if 'slt' in b else np.zeros(len(idx))
    print(f'  {name:16s} {m.sum():6d} {(a["cl"][idx] >= 0.5).sum():11d} '
          f'{int(((d >= FLAKE_D_MIN) & (d <= FLAKE_D_MAX)).sum()):11d} '
          f'{np.median(slt_out):9.1f}')
rest = flip[~claimed]
if rest.size:
    d = a['dl'][rest]
    print(f'  {"elsewhere":16s} {len(rest):6d} {(a["cl"][rest] >= 0.5).sum():11d} '
          f'{int(((d >= FLAKE_D_MIN) & (d <= FLAKE_D_MAX)).sum()):11d} '
          f'{np.median(b["slt"][rest]):9.1f}')

# ------------------------------------------------------------ 4. the recommendation
print('\n' + '=' * 100)
print('\nWHAT A LAKE FILL WOULD HAVE TO SUPPLY\n')
need = flip[a['cl'][flip] < 0.5]
print(f'  {len(flip) - len(need)} cells already carry cl >= 0.5 in the input: '
      f'setting cl = 1 there\n  restores what the input said.')
print(f'  {len(need)} cells were open ocean with no lake cover: a lake fill has to\n'
      f'  INVENT both cl and dl for them.')
if len(need):
    d = a['dl'][need]
    print(f'    their pristine dl: min {d.min():.2f} median {np.median(d):.2f} '
          f'max {d.max():.2f}')
    print(f'    inside FLake\'s band as-is: '
          f'{int(((d >= FLAKE_D_MIN) & (d <= FLAKE_D_MAX)).sum())} of {len(need)}')
print()
