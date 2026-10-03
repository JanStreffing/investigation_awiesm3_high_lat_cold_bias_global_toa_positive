"""Which sea-ice knobs are hemispherically selective BY PHYSICS, not by fiat?

Hemispheric tuning parameters are forbidden in this campaign, and rightly. But the two
polar packs are physically different surfaces, so some knobs act almost entirely in one
hemisphere without anyone having to say so:

  * MELT PONDS.  ice_meltponds.F90:125 disables ponds wherever snow exceeds hs1 = 3 cm
    (also below hi_min ice).  Antarctic ice carries snow through its melt season, Arctic
    ice does not, so rfracmax / pndaspect / albpnd are close to NH-only.
  * SNOW vs BARE-ICE ALBEDO.  ice_thermo_cpl.F90 weights albsn/albsnm against albi/albim
    by snow cover.  If the SH melt-season pack is snow-covered and the NH is not, then
    albsn/albsnm are an SH lever and albi/albim an NH lever, from the same namelist.

This measures both, over the melt season of each hemisphere, on the pack only
(a_ice >= CONC), so the split can be checked instead of assumed.

Usage:  Y0=2030 Y1=2039 python3 scripts/analysis/ice_hemispheric_levers.py
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PICAL_momixoff')
Y0, Y1 = int(os.environ.get('Y0', 2030)), int(os.environ.get('Y1', 2039))
CONC = float(os.environ.get('CONC', 0.5))
P = f'{R}/{ARM}/outdata/fesom'
VARS = ('a_ice', 'h_snow', 'apnd', 'hpnd', 'h_ice')


# Some fields are written daily (365/366 records) and some monthly (12), and leap years
# change the daily length, so everything is reduced to 12 monthly means before stacking.
_ML = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


def _to_monthly(a):
    if a.shape[0] == 12:
        return a
    leap = a.shape[0] == 366
    ml = list(_ML)
    if leap:
        ml[1] = 29
    out, i = [], 0
    for n in ml:
        out.append(np.nanmean(a[i:i + n], 0)); i += n
    return np.stack(out)


def load(v):
    acc = []
    for y in range(Y0, Y1 + 1):
        try:
            with xr.open_dataset(f'{P}/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
                acc.append(_to_monthly(np.squeeze(d[v].values)))
        except Exception:
            pass
    return (np.nanmean(np.stack(acc), 0), ) if acc else (None, )


M = {v: load(v)[0] for v in VARS}
with xr.open_dataset(f'{P}/a_ice.fesom.gr.{Y0}.nc', decode_times=False) as d:
    lat = np.squeeze(d['lat'].values)
LAT = np.broadcast_to(lat[:, None], M['a_ice'].shape[1:])
W = np.cos(np.deg2rad(LAT))

# melt season of each hemisphere: NH Jun-Aug, SH Dec-Feb
SEAS = {'NH melt (JJA), 55-90N': ([5, 6, 7], 55, 90),
        'SH melt (DJF), 90-55S': ([11, 0, 1], -90, -55)}

print(__doc__.split('Usage:')[0])
print(f'{ARM}  {Y0}-{Y1}, pack only (a_ice >= {CONC})\n')
print(f'{"":<24}{"h_snow":>9}{"snow-free":>11}{"apnd":>8}{"hpnd":>8}{"h_ice":>8}')
print(f'{"":<24}{"[m]":>9}{"% of pack":>11}{"frac":>8}{"[m]":>8}{"[m]":>8}')
for name, (mon, la, lb) in SEAS.items():
    ai = np.nanmean(np.stack([M['a_ice'][m] for m in mon]), 0)
    k = (ai >= CONC) & (LAT >= la) & (LAT <= lb) & np.isfinite(ai)
    if k.sum() < 10:
        print(f'  {name:<22} fewer than 10 cells'); continue
    av = lambda f: float(np.average(np.nan_to_num(np.nanmean(np.stack([f[m] for m in mon]), 0))[k],
                                    weights=W[k]))
    hs = np.nan_to_num(np.nanmean(np.stack([M['h_snow'][m] for m in mon]), 0))
    snowfree = float(np.average((hs[k] < 0.03).astype(float), weights=W[k])) * 100
    print(f'  {name:<22}{av(M["h_snow"]):9.3f}{snowfree:11.1f}{av(M["apnd"]):8.3f}'
          f'{av(M["hpnd"]):8.3f}{av(M["h_ice"]):8.2f}')

print("""
  snow-free % is the fraction of the pack below the 3 cm pond threshold hs1, i.e. where
  ponds are allowed to exist at all and where albi/albim rather than albsn/albsnm set the
  pond-free albedo.  A large NH / small SH split means the namelist already separates the
  hemispheres through the physics, and no hemispheric parameter is needed.
""")
