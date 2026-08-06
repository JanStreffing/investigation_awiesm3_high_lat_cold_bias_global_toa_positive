"""The freeze-up latent-heat buffer: is the winter soil collapse a MOISTURE problem?

WHY THIS EXISTS.  The snow-cover schemes (N2/O1/O2) cost 17-25 K of Siberian winter
soil temperature while leaving DJF snow cover unchanged to 0.002 and DJF SWE
unchanged to 1 kg/m2 (snow_state_diag.py).  Cover cannot be the winter channel.
The competing explanation is thermal INERTIA rather than thermal resistance.

Soil water freezing is the dominant winter heat buffer in permafrost soils.  For the
top layer alone (D1 = 0.07 m) at theta = 0.3 m3/m3 the frozen-water latent heat is
    0.07 * 0.3 * 1000 * 334 kJ = 7.0 MJ/m2
against a dry-soil heat capacity of 2.5 MJ/m3/K * 0.07 = 0.175 MJ/m2/K, i.e. the
latent heat is worth ~40 K of layer-1 temperature.  A scheme that dries the soil
before freeze-up removes that buffer and the soil free-falls to air temperature --
which is exactly what N2 does (soil 240.4 K against a 238.6 K snowpack, a 2 K
gradient under half a metre of snow, physically impossible by conduction alone).

The schemes plausibly dry the soil: they change melt timing and the partition of
melt water between infiltration and runoff, and they change how much of the box is
snow tile (which evaporates as sublimation) versus bare/vegetated tile (which
transpires).  None of that shows up in SWE or cover.

WHAT TO READ
  * swvl1..4 in Sep/Oct -- the state at freeze-up.  If N2/O1/O2 are drier and the
    drying ORDERS the runs the way the soil damage does (N2 worst, O2 least), the
    buffer is the mechanism.
  * the latent-heat buffer in MJ/m2 and its K-equivalent, so the size of the effect
    can be compared with the 17-25 K that has to be explained.
  * evaporation (e) and runoff (ro) as the two routes by which the water could have
    left, to say WHICH one the scheme opened.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

RUNS = [('N1 ref off', 'amip_N1_snowdiag'),
        ('N2 cur',     'amip_N2_snowdiag_scf'),
        ('O1 m4',      'amip_O1_scf_m4'),
        ('O2 m3',      'amip_O2_scf_m3')]
YEARS = range(1876, 1916)
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()

DZ = [0.07, 0.21, 0.72, 1.89]          # HTESSEL layer thicknesses [m]
RHOW, LFUS = 1000.0, 3.34e5            # [kg/m3], [J/kg]
CSOIL = 2.5e6                          # dry soil volumetric heat capacity [J/m3/K]


def monthly(run, var, sel):
    """(12, ncell) climatology from monthly files."""
    acc, n = None, 0
    for y in YEARS:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values
        if a.shape[0] != 12:
            continue
        a = a.reshape(12, -1)[:, sel]
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, n) if n else (None, 0)


with xr.open_dataset(LSMF) as d:
    lsm = d['lsm'].isel(time_counter=0).values
    lat, lon = d['lat'].values, d['lon'].values
LA = np.broadcast_to(lat[:, None], lsm.shape)
LO = np.broadcast_to(lon[None, :], lsm.shape)
boxm = (LA >= 55) & (LA <= 75) & (LO >= 60) & (LO <= 180) & (lsm > 0.5)
sel = np.flatnonzero(boxm.ravel())
w = np.cos(np.deg2rad(LA)).ravel()[sel]
W = w / w.sum()

print(__doc__.split('WHAT TO READ')[0])
print('=' * 92)
print(f'Siberian land 55-75N 60-180E, {sel.size} cells\n')

V = {}
for lab, run in RUNS:
    d = {}
    for v in ('swvl1', 'swvl2', 'swvl3', 'swvl4', 'stl1', 'e', 'ro'):
        a, n = monthly(run, v, sel)
        if a is not None:
            d[v] = a @ W
    if 'swvl1' not in d:
        print(f'  {lab}: no monthly output'); continue
    V[lab] = d
    print(f'  {lab:12s} {n} yr')

ref = V['N1 ref off']

for v in ('swvl1', 'swvl2', 'swvl3', 'swvl4'):
    print(f'\n\n{v.upper()} -- volumetric soil water [m3/m3]\n')
    print(f'  {"run":12s}' + ''.join(f'{m:>8s}' for m in MON))
    for lab in V:
        print(f'  {lab:12s}' + ''.join(f'{x:8.4f}' for x in V[lab][v]))
    print(f'  {"-- delta":12s}')
    for lab in [l for l in V if l != 'N1 ref off']:
        print(f'  {lab:12s}' + ''.join(f'{x:+8.4f}' for x in V[lab][v] - ref[v]))

# ---- the buffer, in the units that matter -------------------------------------
print('\n\n  LATENT-HEAT BUFFER at freeze-up (Sep-Oct mean), all four layers\n')
print(f'  {"run":12s}{"water [kg/m2]":>15s}{"L [MJ/m2]":>12s}{"K-equiv":>10s}{"dK vs N1":>10s}')
def buffer(d):
    m = sum(0.5 * (d[f'swvl{i+1}'][8] + d[f'swvl{i+1}'][9]) * DZ[i] * RHOW for i in range(4))
    return m, m * LFUS / 1e6, m * LFUS / (CSOIL * sum(DZ))
b0 = buffer(ref)
for lab in V:
    m, mj, kq = buffer(V[lab])
    print(f'  {lab:12s}{m:15.1f}{mj:12.2f}{kq:10.1f}{kq - b0[2]:+10.2f}')

print('\n\n  WHERE THE WATER WENT -- evaporation and runoff, JJA+SON totals\n')
print(f'  {"run":12s}{"e JJA":>10s}{"e SON":>10s}{"ro JJA":>10s}{"ro SON":>10s}')
for lab in V:
    if 'e' not in V[lab]:
        continue
    f = lambda v, ms: sum(V[lab][v][m] for m in ms)
    print(f'  {lab:12s}{f("e",[5,6,7]):10.4f}{f("e",[8,9,10]):10.4f}'
          f'{f("ro",[5,6,7]):10.4f}{f("ro",[8,9,10]):10.4f}')

print("""
  READING IT.  The K-equivalent column is the whole-column temperature change the
  frozen water can absorb before the soil is free to cool further.  If dK vs N1 is
  a fraction of a kelvin, moisture is NOT the mechanism and this hypothesis dies
  with the others.  If it is several kelvin AND ordered N2 < O1 < O2, it is.""")
