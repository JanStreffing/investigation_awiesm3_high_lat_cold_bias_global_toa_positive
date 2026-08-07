"""Round 20 verdict: does the fitted depletion cure the winter it was built to cure?

THE PRE-REGISTERED FALSIFIERS, written before P3/P4 ran:
  1. DJF soil (stl1/stl2) must return to the N1/K1 reference.  If it does not, the
     distribution diagnosis is wrong and the whole round-20 argument fails.
  2. January f_full -- the area fraction at cover >= 0.999 -- must return to ~0.96
     from N2's 0.773.  That is the quantity that actually predicted the damage
     (per-cell correlation +0.75 against +0.32 for mean cover).
  3. For P4 specifically: if it recovers the JJA gain BUT reopens the winter soil
     bias, then spring depletion and winter insulation really are coupled at grid
     scale, and the flat frontier found on station data does not survive.

WHAT EACH RUN IS
  N1  K1 base, scheme OFF          -- the reference, and the run that is RIGHT vs
                                      RIHMI station soil (+0.8 K at 0.2 m)
  N2  tanh, as-tuned               -- the broken one, -20.2 K vs the same stations
  P3  fitted curve, SCALE=1        -- observations taken literally
  P4  fitted curve, SCALE=3        -- d_c x3 for 100 km sub-grid variance,
                                      the UNCALIBRATED bracket

Cover is reconstructed per day per cell from each run's OWN formula and only then
averaged -- averaging depth first and applying the nonlinear curve after is the
error that voided snow_daily_diag.py, and it is not repeated here.

Mode 3 needs the high-vegetation fraction cvh, which the N/O diagnostics never
loaded; the P runs carry it as daily output.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

OBSD = '/work/ab0246/a270092/obs/RIHMI-WDC/data'
YEARS = range(1876, 1916)
DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
DOY_MONTH = np.repeat(np.arange(12), DPM)

# mode 3 parameters, mirroring surfbc_ctl_mod.F90 INCLUDING the SWEMIN mass floor.
M3 = dict(dcl=0.014, dch=0.026, md=4.70, bl=1.46, bh=0.40,
          dcmax=0.30, rhoref=200.0, swemin=3.0)
# P5/P6 ran on the PRE-DCMAX-FIX binary, where MIN(DCMAX, MAX(fitted, floor)) let the
# 0.30 m cap clip the floor.  Reconstructing them with the CURRENT (fixed) formula
# would not describe what actually ran, so the old order is reproduced here.
RUNS = [('N1 off',  'amip_N1_snowdiag',     ('linear', None)),
        ('N2 tanh', 'amip_N2_snowdiag_scf', ('tanh',   (0.016, 100.0, 1.6))),
        ('P3 sw3',  'amip_P3_scffit',       ('mode3',  dict(M3, scale=1.0, swemin=3.0))),
        ('P5 sw15', 'amip_P5_swemin15',     ('mode3',  dict(M3, scale=1.0, swemin=15.0))),
        ('P6 sw30', 'amip_P6_swemin30',     ('mode3',  dict(M3, scale=1.0, swemin=30.0)))]
# Rutgers 24 km, Siberian box (rutgers_vs_model_cover.py); the round-21 target
RUTGERS = {8: 0.067, 9: 0.599, 10: 0.969, 11: 0.997, 0: 0.999, 1: 1.000,
           2: 0.999, 3: 0.956, 4: 0.655, 5: 0.202}


def scf(kind, p, depth, rho, cvh):
    if kind == 'linear':
        return np.clip(10.0 * depth, 0.0, 1.0)
    live = (depth > 1e-6) & (depth * rho > 1e-6)
    if kind == 'tanh':
        z0, rn, m = p
        s = 2.5 * z0 * (np.maximum(rho, 50.0) / rn) ** m
        return np.where(live, np.clip(np.tanh(depth / np.maximum(s, 1e-6)), 0, 1), 0.0)
    r = np.maximum(rho, 50.0) / p['rhoref']
    floor = p['swemin'] / np.maximum(rho, 1.0)
    # AS RUN: DCMAX caps the floor too (the bug fixed in 562df81, after these runs).
    dcl = np.minimum(p['dcmax'], np.maximum(p['scale'] * p['dcl'] * r ** p['md'], floor))
    dch = np.minimum(p['dcmax'], np.maximum(p['scale'] * p['dch'] * r ** p['md'], floor))
    sl = np.clip((depth / np.maximum(dcl, 1e-9)) ** p['bl'], 0, 1)
    sh = np.clip((depth / np.maximum(dch, 1e-9)) ** p['bh'], 0, 1)
    f = np.clip(cvh, 0.0, 1.0)
    return np.where(live, np.clip((1 - f) * sl + f * sh, 0, 1), 0.0)


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
box = (LA >= 55) & (LA <= 75) & (LO >= 60) & (LO <= 180) & (lsm > 0.5)
sel = np.flatnonzero(box.ravel())
w = np.cos(np.deg2rad(LA)).ravel()[sel]; W = w / w.sum()

print(__doc__)
print('=' * 100)
print(f'Siberian land 55-75N 60-180E, {sel.size} cells, {len(list(YEARS))} yr target\n')

S = {}
for lab, run, (kind, p) in RUNS:
    acc = {k: np.zeros(12) for k in ('cover', 'full', 'swe', 'depth', 'rho', 'stl1', 'stl2', 'tsn')}
    n = 0
    for y in YEARS:
        sd, rsn = load(run, 'sd', y, sel), load(run, 'rsn', y, sel)
        if sd is None or rsn is None:
            continue
        cvh = load(run, 'cvh', y, sel)
        if cvh is None:
            cvh = np.zeros_like(sd)
        rho = np.maximum(rsn, 1e-6)
        swe = sd * 1000.0
        depth = swe / rho
        c = scf(kind, p, depth, rho, cvh)
        st1, st2, tsn = (load(run, v, y, sel) for v in ('stl1', 'stl2', 'tsn'))
        for m in range(12):
            k = DOY_MONTH == m
            acc['cover'][m] += (c[k].mean(axis=0) @ W)
            acc['full'][m] += ((c[k] >= 0.999).mean(axis=0) @ W)
            acc['swe'][m] += (swe[k].mean(axis=0) @ W)
            acc['depth'][m] += (depth[k].mean(axis=0) @ W)
            acc['rho'][m] += (rho[k].mean(axis=0) @ W)
            for v, a in (('stl1', st1), ('stl2', st2), ('tsn', tsn)):
                if a is not None:
                    acc[v][m] += (a[k].mean(axis=0) @ W)
        n += 1
    if not n:
        print(f'  {lab}: NO DAILY OUTPUT'); continue
    S[lab] = {k: v / n for k, v in acc.items()}; S[lab]['n'] = n
    print(f'  {lab:9s} {n} yr')

ref = S['N1 off']
DJF = [11, 0, 1]; JJA = [5, 6, 7]
mn = lambda a, ms: float(np.mean([a[m] for m in ms]))

# ---------------------------------------------------------------- falsifier 2 --
print('\n\n1. f_full -- area fraction at cover >= 0.999   [FALSIFIER 2: P3/P4 Jan ~0.96]\n')
print(f'  {"run":9s}' + ''.join(f'{MON[m]:>8s}' for m in (8, 9, 10, 11, 0, 1, 2, 3, 4)))
for lab in S:
    print(f'  {lab:9s}' + ''.join(f'{S[lab]["full"][m]:8.3f}' for m in (8, 9, 10, 11, 0, 1, 2, 3, 4)))
print(f'\n  box-mean cover:')
for lab in S:
    print(f'  {lab:9s}' + ''.join(f'{S[lab]["cover"][m]:8.3f}' for m in (8, 9, 10, 11, 0, 1, 2, 3, 4)))

# ---------------------------------------------------------------- falsifier 1 --
print('\n\n2. SOIL TEMPERATURE, delta vs N1 [K]   [FALSIFIER 1: must return to ~0]\n')
for v in ('stl1', 'stl2'):
    print(f'  --- {v}')
    print(f'  {"run":9s}' + ''.join(f'{MON[m]:>8s}' for m in (9, 10, 11, 0, 1, 2)) + f'{"DJF":>9s}')
    print(f'  {"N1 (abs)":9s}' + ''.join(f'{ref[v][m]-273.15:8.1f}' for m in (9, 10, 11, 0, 1, 2))
          + f'{mn(ref[v],DJF)-273.15:9.1f}')
    for lab in [l for l in S if l != 'N1 off']:
        d = [S[lab][v][m] - ref[v][m] for m in (9, 10, 11, 0, 1, 2)]
        print(f'  {lab:9s}' + ''.join(f'{x:+8.2f}' for x in d)
              + f'{mn(S[lab][v],DJF)-mn(ref[v],DJF):+9.2f}')

# ------------------------------------------------------------------ the spring --
print('\n\n2b. COVER vs RUTGERS 24 km  [round-21 target: cut the Sep/Oct excess]\n')
print(f'  {"run":9s}' + ''.join(f'{MON[m]:>8s}' for m in (8, 9, 10, 11, 0, 3, 4)))
print(f'  {"Rutgers":9s}' + ''.join(f'{RUTGERS[m]:8.3f}' for m in (8, 9, 10, 11, 0, 3, 4)))
for lab in S:
    print(f'  {lab:9s}' + ''.join(f'{S[lab]["cover"][m]-RUTGERS[m]:+8.3f}' for m in (8, 9, 10, 11, 0, 3, 4)))
print('  (rows after Rutgers are model MINUS observed)')

print('\n\n3. SPRING DEPLETION -- does the fitted curve buy any?   [what the scheme is FOR]\n')
print(f'  {"run":9s}{"Apr cov":>9s}{"May cov":>9s}{"Jun cov":>9s}{"May SWE":>9s}{"Jun SWE":>9s}')
for lab in S:
    print(f'  {lab:9s}' + ''.join(f'{S[lab]["cover"][m]:9.3f}' for m in (3, 4, 5))
          + ''.join(f'{S[lab]["swe"][m]:9.2f}' for m in (4, 5)))

# ----------------------------------------------------------- the pack and gradient
print('\n\n4. PACK STATE and SOIL-TO-SNOW GRADIENT (Jan)   [N1 22 K, N2 2 K -- unphysical]\n')
print(f'  {"run":9s}{"SWE":>9s}{"depth":>9s}{"rho":>9s}{"tsn":>9s}{"stl1":>9s}{"grad":>9s}')
for lab in S:
    g = S[lab]['stl1'][0] - S[lab]['tsn'][0]
    print(f'  {lab:9s}{S[lab]["swe"][0]:9.1f}{S[lab]["depth"][0]:9.3f}{S[lab]["rho"][0]:9.1f}'
          f'{S[lab]["tsn"][0]:9.1f}{S[lab]["stl1"][0]:9.1f}{g:9.1f}')

# --------------------------------------------------- observational check vs RIHMI
with xr.open_dataset(f'{OBSD}/RIHMI-WDC_tpg.nc') as d:
    la, lo = d['lat'].values, d['lon'].values
    k = (la >= 55) & (la <= 75) & (lo >= 60) & (lo <= 180)
    ts = d['tsoil'].isel(station=k).values
    dep = d['depth'].values
    tm = d['time'].values
mo = tm.astype('datetime64[M]').astype(int) % 12
i20 = int(np.abs(dep - 0.2).argmin())
obs = np.nanmean([np.nanmean(ts[:, mo == m, i20]) for m in DJF])
print(f'\n\n5. vs RIHMI STATION SOIL at {dep[i20]:.2f} m, DJF   [N1 is RIGHT to +0.8 K]\n')
print(f'  observed (1963-2024, 43 stations): {obs:+6.1f} degC')
print(f'  {"run":9s}{"stl2 DJF":>11s}{"bias":>8s}')
for lab in S:
    v = mn(S[lab]['stl2'], DJF) - 273.15
    print(f'  {lab:9s}{v:11.1f}{v-obs:+8.1f}')
print('  (model is PI, obs 1963-2024; Siberian winter warming makes the PI-equivalent')
print('   observation ~3 K colder, so read N1 as ~+3 K warm and the tanh as ~-17 K cold.)')

print("""

  VERDICT RULES.  P3 passes falsifier 1 if its DJF soil delta is within the DJF
  detection threshold (+-0.588 K) of N1, and falsifier 2 if January f_full is back
  near 0.96.  If P3 passes both but shows no spring depletion (section 3), the
  scheme is SAFE but BUYS NOTHING at SCALE=1, and everything then rests on P4 --
  whose SCALE is uncalibrated and must not be adopted on this evidence alone.""")
