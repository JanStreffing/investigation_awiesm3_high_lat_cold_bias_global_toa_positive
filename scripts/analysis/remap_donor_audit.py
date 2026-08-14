"""Which soil does each remapped LPJ-GUESS cell inherit from, and does it fit?

WHY THIS EXISTS.  state_remap gives a target gridcell the state of its nearest
source cell.  The state carries water, ice and below-wilting-point water as
absolute volumetric fractions; it does NOT carry porosity, awc or Fpwp_ref,
which are recomputed from the TARGET's soil code.  So a donor on a different
soil overfills its new home, and the run dies on day 0 -- first on Frac_air,
then on wcont, then on the freeze/thaw mass balance.  Clamping each check in
turn is whack-a-mole; the question is whether the donated state is compatible
at all.

THE MEASUREMENT.  For every cell the TCO95-land gridlist adds over L096, find
its nearest source cell the same way guessserializer.cpp does (flat spherical
approximation, cos-weighted longitude) and compare soil codes.  That says:

  * how many new cells there are and how far their donors sit,
  * how many inherit from a DIFFERENT soil code, which is the failure set,
  * whether a same-soil donor exists nearby, which decides whether restricting
    the donor search is a viable fix or whether the state must be respun.

A TRAP THIS AVOIDS.  TCO95-land.msk is stored INVERTED (ocean = 1), as
eceframework.cpp says; L096.msk is not.  Getting that backwards silently swaps
source and target and makes the whole audit meaningless.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, netCDF4, xarray as xr

OCP = '/work/ab0246/a270092/software/ocp-tool/output/TCO95_CORE3'
MASKS = f'{OCP}/oasis_mct3_input/masks.nc'
GRIDS = f'{OCP}/oasis_mct3_input/grids.nc'
SLT_OLD = '/work/ab0246/a270092/input/lpj-guess/slt/slt_TCO95_CORE3.nc'      # 547
SLT_NEW = '/work/ab0246/a270092/input/lpj-guess/slt/slt_TCO95_CORE3_v2.nc'   # 233

print(__doc__)
print('=' * 96)

with netCDF4.Dataset(MASKS) as d:
    L = np.asarray(d.variables['L096.msk'][:]).ravel()          # 1 = land
    T = np.asarray(d.variables['TCO95-land.msk'][:]).ravel()    # INVERTED: 1 = ocean
with netCDF4.Dataset(GRIDS) as d:
    lon = np.asarray(d.variables['A096.lon'][:]).ravel()
    lat = np.asarray(d.variables['A096.lat'][:]).ravel()

src = L == 1          # the gridlist the spin-up state was built on
tgt = T == 0          # the gridlist the mask fix asks for
new = tgt & ~src      # cells with no state of their own
print(f'\n  source (L096)      : {src.sum()}')
print(f'  target (TCO95-land): {tgt.sum()}')
print(f'  NEW cells          : {new.sum()}      dropped: {int((src & ~tgt).sum())}')

for tag, path in (('old 547', SLT_OLD), ('v2  233', SLT_NEW)):
    ds = xr.open_dataset(path)
    slt = np.squeeze(ds[[k for k in ds.data_vars][0]].values).astype(int)

    si = np.where(src)[0]
    ni = np.where(new)[0]
    slon, slat = lon[si], lat[si]

    donor = np.empty(len(ni), dtype=int)
    dist = np.empty(len(ni))
    for k, i in enumerate(ni):
        dlon = (lon[i] - slon) * np.cos(np.deg2rad((lat[i] + slat) * 0.5))
        dlat = lat[i] - slat
        d2 = dlon * dlon + dlat * dlat
        j = int(np.argmin(d2))
        donor[k] = si[j]
        dist[k] = np.sqrt(d2[j])

    t_slt, d_slt = slt[ni], slt[donor]
    same = t_slt == d_slt
    print(f'\n  --- soil map {tag} ---')
    print(f'  NN distance  : median {np.median(dist):.3f} deg  max {dist.max():.3f} deg')
    print(f'  donor soil == target soil : {int(same.sum())} of {len(ni)} '
          f'({100*same.mean():.0f} %)')
    print(f'  MISMATCHED                : {int((~same).sum())}   <- the failure set')
    if (~same).any():
        import collections
        c = collections.Counter(zip(d_slt[~same].tolist(), t_slt[~same].tolist()))
        print(f'  top donor->target transitions:')
        for (a, b), n in c.most_common(6):
            print(f'      soil {a} -> {b}   {n} cells')
    # is a same-soil donor available nearby?
    if (~same).any():
        worst = []
        for k in np.where(~same)[0]:
            i = ni[k]
            cand = si[slt[si] == t_slt[k]]
            if cand.size == 0:
                worst.append(np.inf); continue
            dlon = (lon[i] - lon[cand]) * np.cos(np.deg2rad((lat[i] + lat[cand]) * 0.5))
            dlat = lat[i] - lat[cand]
            worst.append(np.sqrt(np.min(dlon * dlon + dlat * dlat)))
        worst = np.array(worst)
        fin = np.isfinite(worst)
        print(f'  nearest SAME-soil donor   : median {np.median(worst[fin]):.2f} deg  '
              f'max {worst[fin].max():.2f} deg  '
              f'({int((~fin).sum())} cells with no same-soil donor anywhere)')
print()
