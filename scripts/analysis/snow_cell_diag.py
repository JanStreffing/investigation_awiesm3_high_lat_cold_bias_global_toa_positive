"""PER-CELL snow cover vs soil damage: does the box mean hide the mechanism?

WHY THIS EXISTS.  snow_state_diag.py established two facts that no box-mean story
survives:
  * DJF cover is the SAME in all four runs to within 0.002, and the soil differs
    by up to 23 K.
  * October soil damage is -4.94 / -5.02 / -4.95 K in N2 / O1 / O2 while October
    cover differs by 0.10 between them.  Cover moves, damage does not.
and one that inverts the physics as reconstructed on box means:
  * N1 carries a 22 K gradient from soil (265.0 K) to snow (242.7 K); N2 carries
    2 K (240.4 / 238.6).  The soil is far MORE tightly tied to the pack in N2,
    while the box-mean conductance says it should be LESS.

The remaining explanation is spatial.  tanh(d/L) and min(1,10d) can deliver the
same box-mean cover from completely different per-cell distributions -- the linear
ramp CLIPS, so every cell with more than 10 cm reports exactly 1.0 and the box mean
falls short of 1 only because some cells are bare; the tanh never reports 1.0
anywhere, so the same box mean is built from EVERY cell being slightly uncovered.
Those are physically opposite states and the soil response to cover is strongly
nonlinear (a cell that goes bare in polar night crashes; a cell at 0.99 does not).

This script therefore stops averaging over the box before it has looked at it.

WHAT TO READ
  * cover HISTOGRAM per run per month: is N1's cover a spike at 1.0 with a bare
    tail, while N2/O1/O2 are a broad hump just below 1?  That is the whole thesis.
  * f_full: area fraction of the box with cover >= 0.999.  If N1 is ~0.95 and the
    scheme runs are ~0.00, the schemes have removed complete snow cover from
    EVERY cell without changing the mean, and that is the mechanism.
  * per-cell correlation of dCover with dSoil, and dSoil binned by dCover: if the
    damage is concentrated where cover fell, cover is the channel after all and
    only the box mean was lying.  If dSoil is uniform and uncorrelated with
    dCover, it is not.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

RUNS = [('N1 ref off', 'amip_N1_snowdiag',     None),
        ('N2 cur',     'amip_N2_snowdiag_scf', (0.016, 100.0, 1.6)),
        ('O1 m4',      'amip_O1_scf_m4',       (0.018, 170.0, 4.0)),
        ('O2 m3',      'amip_O2_scf_m3',       (0.018, 170.0, 3.0))]
YEARS = range(1876, 1896)          # 20 yr is ample for a spatial-structure question
DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
DOY_MONTH = np.repeat(np.arange(12), DPM)
SHOW = [8, 9, 10, 11, 0, 1]        # Sep..Feb -- where the damage is built
RQSNCR_INV = 10.0


def cover(depth, rho, p):
    if p is None:
        return np.clip(RQSNCR_INV * depth, 0.0, 1.0)
    z0, rho_new, m = p
    scale = 2.5 * z0 * (np.maximum(rho, 50.0) / rho_new) ** m
    c = np.tanh(depth / np.maximum(scale, 1e-6))
    return np.where((depth > 1e-6) & (depth * rho > 1e-6), np.clip(c, 0, 1), 0.0)


def load(run, var, y, sel):
    for pat in (f'atm_remapped_1d_{var}_1d_{y}-{y}.nc', f'atm_remapped_1d_{var}_{y}-{y}.nc'):
        f = f'{RT}/{run}/outdata/oifs/{pat}'
        if os.path.exists(f):
            with xr.open_dataset(f, decode_times=False) as d:
                a = d[var].values
            if a.shape[0] == 366:
                a = np.delete(a, 59, axis=0)
            return a.reshape(a.shape[0], -1)[:, sel] if a.shape[0] == 365 else None
    return None


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
print('=' * 96)
print(f'Siberian land 55-75N 60-180E, {sel.size} cells, {len(list(YEARS))} yr\n')

# per-cell monthly means of cover / full-cover fraction / stl1
C = {}
for lab, run, p in RUNS:
    cov = np.zeros((12, sel.size)); full = np.zeros((12, sel.size))
    soil = np.zeros((12, sel.size)); n = 0
    for y in YEARS:
        sd, rsn, st = (load(run, v, y, sel) for v in ('sd', 'rsn', 'stl1'))
        if sd is None or rsn is None or st is None:
            continue
        rho = np.maximum(rsn, 1e-6)
        c = cover(sd * 1000.0 / rho, rho, p)
        for m in range(12):
            k = DOY_MONTH == m
            cov[m] += c[k].mean(axis=0)
            full[m] += (c[k] >= 0.999).mean(axis=0)
            soil[m] += st[k].mean(axis=0)
        n += 1
    if not n:
        print(f'  {lab}: no daily output'); continue
    C[lab] = dict(cov=cov / n, full=full / n, soil=soil / n, n=n)
    print(f'  {lab:12s} {n} yr')

ref = C['N1 ref off']

# ---- 1. the decisive one: how much of the box is at COMPLETE cover? -----------
print('\n\n1. f_full -- area fraction of the box with cover >= 0.999 (time+area weighted)\n')
print(f'  {"run":12s}' + ''.join(f'{MON[m]:>9s}' for m in SHOW))
for lab in C:
    print(f'  {lab:12s}' + ''.join(f'{C[lab]["full"][m] @ W:9.3f}' for m in SHOW))
print('\n  Box-mean cover, for comparison:')
for lab in C:
    print(f'  {lab:12s}' + ''.join(f'{C[lab]["cov"][m] @ W:9.3f}' for m in SHOW))

# ---- 2. the cover histogram ---------------------------------------------------
EDGES = [0.0, 0.1, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999, 1.001]
for m in (9, 11, 0):
    print(f'\n\n2. COVER HISTOGRAM, {MON[m]} -- area fraction of cells in each bin\n')
    hdr = ''.join(f'{f"{EDGES[i]:g}-{EDGES[i+1]:g}":>10s}' for i in range(len(EDGES) - 1))
    print(f'  {"run":12s}{hdr}')
    for lab in C:
        c = C[lab]['cov'][m]
        row = [W[(c >= EDGES[i]) & (c < EDGES[i + 1])].sum() for i in range(len(EDGES) - 1)]
        print(f'  {lab:12s}' + ''.join(f'{v:10.3f}' for v in row))

# ---- 3. is the damage where the cover fell? ----------------------------------
print('\n\n3. PER-CELL dSoil vs dCover  (Sep-Nov cover change -> DJF soil change)\n')
autumn = lambda a: (a[8] + a[9] + a[10]) / 3.0
djf = lambda a: (a[0] + a[1] + a[11]) / 3.0
print(f'  {"run":12s}{"corr":>8s}{"dSoil<->dFull":>15s}   dSoil binned by dCover decile (cold->warm)')
for lab in [l for l in C if l != 'N1 ref off']:
    dc = autumn(C[lab]['cov']) - autumn(ref['cov'])
    df = autumn(C[lab]['full']) - autumn(ref['full'])
    ds = djf(C[lab]['soil']) - djf(ref['soil'])
    r1 = np.corrcoef(dc, ds)[0, 1]
    r2 = np.corrcoef(df, ds)[0, 1]
    q = np.quantile(dc, np.linspace(0, 1, 6))
    binned = [ds[(dc >= q[i]) & (dc <= q[i + 1])].mean() for i in range(5)]
    print(f'  {lab:12s}{r1:+8.2f}{r2:+15.2f}   ' + ' '.join(f'{v:+7.1f}' for v in binned))

# ---- 4. is the damage uniform or concentrated? -------------------------------
print('\n\n4. DISTRIBUTION of per-cell DJF dSoil [K] -- uniform shift or a few crashed cells?\n')
print(f'  {"run":12s}' + ''.join(f'{f"p{p}":>9s}' for p in (1, 5, 25, 50, 75, 95, 99)) + f'{"mean":>9s}')
for lab in [l for l in C if l != 'N1 ref off']:
    ds = djf(C[lab]['soil']) - djf(ref['soil'])
    print(f'  {lab:12s}' + ''.join(f'{np.percentile(ds, p):9.1f}' for p in (1, 5, 25, 50, 75, 95, 99))
          + f'{ds @ W:9.1f}')

print("""
  READING IT.  If f_full collapses from ~1 in N1 to ~0 in the scheme runs while the
  MEAN cover is unchanged, then the schemes have taken every cell off complete cover
  -- a state the as-released model never produces -- and the soil damage follows the
  loss of complete cover, not the loss of mean cover.  That is a structural defect of
  using an asymptotic function for a fraction that must be able to reach 1, and it is
  fixed by a saturation cut (surfbc_ctl_mod.F90:358 does exactly this, commented out),
  not by re-tuning z0/rho_new/m.""")
