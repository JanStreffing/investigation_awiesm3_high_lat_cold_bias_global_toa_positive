"""Daily snow/soil diagnostic: is the winter soil collapse SEEDED IN AUTUMN?

WHY THIS EXISTS. The snow-cover depletion scheme (ECE_SNOW_SCF) buys the largest summer
gain in the campaign and costs a catastrophic winter: coupled DJF soil -16.2 K, January
-18.3 K. Monthly data suggested a mechanism -- the reconstructed cover deficit peaks at
-0.075 in OCTOBER exactly as the soil starts to diverge, then vanishes while the soil runs
on to -24 K, i.e. autumn seeding plus thermal memory, with midwinter cover identical
(-0.0004 in January). But that rested on monthly means, and two things could not be settled:

  1. Does the soil really run away from an October seed, or does something else take over
     in November? Monthly means cannot separate those.
  2. Is the pack RIPE (wet/isothermal) when the spring depletion is needed? Monthly-mean
     tsn is 271.0 K in May while melt-out is late May -- sub-daily melt events are exactly
     what a monthly mean destroys.

Rounds 19 answered by writing DAILY sd/rsn/tsn/asn/stl1/stl2 (per-run file_def override)
for four runs at full campaign length:
    N1  K1 base, scheme OFF          (reference)
    N2  + scheme, current parameters (z0 0.016, rho_new 100, m 1.6)
    O1  + scheme, z0 0.018, rho_new 170, m 4.0
    O2  + scheme, z0 0.018, rho_new 170, m 3.0

O1/O2 re-parameterise the density dependence so that, ON BOX MEANS, October should be
protected (SCF 0.995) while May still depletes (0.495). THAT PREDICTION FAILED: DJF T2m is
-2.633 (O1) and -2.443 (O2) against N2's -2.783 and the K1 reference -0.496. The winter is
not recovered. The pre-registered falsifier fired.

This script asks WHY, at daily resolution, using the one field that was never available
before -- ZCVS reconstructed per day rather than from a monthly mean. Since SCF is strongly
nonlinear in depth, mean(SCF) != SCF(mean), and the daily reconstruction is the honest one.

WHAT TO READ
  * "cover deficit" vs "soil deficit" by pentad through Aug-Dec: if the soil keeps falling
    after the cover difference has closed, seeding+memory survives. If the soil tracks the
    cover day by day, the mechanism is contemporaneous and an autumn-only fix cannot work.
  * "ripe days": fraction of days with tsn >= 273.15 by month. If small in May, a
    melt-state gate would fire too late and lose the spring depletion entirely.
  * O1/O2 vs N2 in the SAME table: did the re-parameterisation actually protect October
    cover as the box-mean arithmetic promised? If October cover is fixed but the soil is
    not, the autumn-seeding hypothesis is wrong, not the parameters.
"""
import numpy as np, xarray as xr, glob, os, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

RUNS = [('N1 ref',   'amip_N1_snowdiag',     None),
        ('N2 cur',   'amip_N2_snowdiag_scf', (0.016, 100.0, 1.6)),
        ('O1 m4',    'amip_O1_scf_m4',       (0.018, 170.0, 4.0)),
        ('O2 m3',    'amip_O2_scf_m3',       (0.018, 170.0, 3.0))]
YEARS = range(1876, 1916)          # skip spin-up; 40 yr of daily data
DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()


def daily(run, var, years=YEARS):
    """Day-of-year climatology (365, ny, nx), leap days dropped."""
    acc, n, lat, lon = None, 0, None, None
    for y in years:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1d_{var}_1d_{y}-{y}.nc'
        if not os.path.exists(f):
            f = f'{RT}/{run}/outdata/oifs/atm_remapped_1d_{var}_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        d = xr.open_dataset(f, decode_times=False)
        a = d[var].values
        lat, lon = d['lat'].values, d['lon'].values
        d.close()
        if a.shape[0] == 366:                      # drop 29 Feb
            a = np.delete(a, 59, axis=0)
        if a.shape[0] != 365:
            continue
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, lat, lon, n) if n else (None, None, None, 0)


def scf(depth, rho, p):
    """The model's own formula, evaluated per day."""
    z0, rho_new, m = p
    L = 2.5 * z0 * (np.maximum(rho, 50.0) / rho_new) ** m
    return np.clip(np.tanh(depth / np.maximum(L, 1e-6)), 0, 1)


lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values
D = {}
for lab, run, p in RUNS:
    sd, lat, lon, n = daily(run, 'sd')
    if sd is None:
        print(f'  {lab}: no daily output'); continue
    rsn, _, _, _ = daily(run, 'rsn')
    tsn, _, _, _ = daily(run, 'tsn')
    st1, _, _, _ = daily(run, 'stl1')
    D[lab] = dict(sd=sd, rsn=rsn, tsn=tsn, stl1=st1, p=p, n=n)

LA = np.broadcast_to(lat[:, None], lsm.shape)
box = (LA >= 55) & (LA <= 75) & \
      np.broadcast_to(((lon >= 60) & (lon <= 180))[None, :], LA.shape) & (lsm > 0.5)
w = np.cos(np.deg2rad(LA))
bm = lambda a: np.average(a[box], weights=w[box])
doy_month = np.repeat(np.arange(12), DPM)

print(__doc__.split('WHAT TO READ')[0])
print('=' * 78)
print(f'Siberian land, day-of-year climatology, {D["N1 ref"]["n"]} yr\n')

# ---- 1. is the pack ripe when the spring depletion is needed? --------------
print('1. RIPE DAYS -- fraction of days with tsn >= 273.15 K (melting point)\n')
print(f'  {"run":8s}' + ''.join(f'{m:>6s}' for m in MON))
for lab in D:
    t = D[lab]['tsn']
    row = []
    for mi in range(12):
        sel = doy_month == mi
        frac = np.average((t[sel][:, box] >= 273.15).mean(axis=0),
                          weights=w[box])
        row.append(frac)
    print(f'  {lab:8s}' + ''.join(f'{v:6.2f}' for v in row))
print('\n  If May is small, a melt-state gate fires LATE and the spring depletion is lost.')

# ---- 2. daily cover deficit vs daily soil deficit --------------------------
print('\n\n2. COVER DEFICIT vs SOIL DEFICIT, by pentad (vs N1) -- the seeding test\n')
ref = D['N1 ref']
depth_ref = ref['sd'] * 1000.0 / np.maximum(ref['rsn'], 1e-6)
old = np.clip(depth_ref * 10.0, 0, 1)                 # as-released linear formula
for lab in [l for l in D if l != 'N1 ref']:
    d = D[lab]
    depth = d['sd'] * 1000.0 / np.maximum(d['rsn'], 1e-6)
    new = scf(depth, d['rsn'], d['p'])
    dcov = np.array([bm(new[i]) - bm(old[i]) for i in range(365)])
    dsoil = np.array([bm(d['stl1'][i]) - bm(ref['stl1'][i]) for i in range(365)])
    print(f'  --- {lab} ---')
    print(f'  {"pentad":>10s}' + ''.join(f'{s:>8s}' for s in
          ('Aug', 'Sep1', 'Sep2', 'Oct1', 'Oct2', 'Nov1', 'Nov2', 'Dec1', 'Dec2', 'Jan')))
    wins = [(212, 243), (243, 258), (258, 273), (273, 288), (288, 303),
            (303, 318), (318, 334), (334, 349), (349, 365), (0, 31)]
    print(f'  {"dCover":>10s}' + ''.join(f'{dcov[a:b].mean():+8.3f}' for a, b in wins))
    print(f'  {"dSoil [K]":>10s}' + ''.join(f'{dsoil[a:b].mean():+8.2f}' for a, b in wins))
    # when does each first exceed a threshold?
    ic = next((i for i in range(212, 365) if dcov[i] < -0.02), None)
    isl = next((i for i in range(212, 365) if dsoil[i] < -1.0), None)
    f = lambda i: '--' if i is None else f'doy {i} ({MON[doy_month[i]]})'
    print(f'  first day cover deficit < -0.02 : {f(ic)}')
    print(f'  first day soil deficit  < -1.0 K: {f(isl)}\n')

print("""  READING IT. If the soil deficit keeps growing through Nov-Jan while the cover
  deficit has already closed, the damage is SEEDED in autumn and carried by thermal
  memory -- an autumn-only fix could work. If the soil tracks the cover day by day, the
  effect is contemporaneous and protecting October alone cannot help.

  And if O1/O2 DID protect October cover (dCover ~ 0 there) but the soil still collapses,
  then the autumn-seeding hypothesis is wrong -- not the parameters.""")
