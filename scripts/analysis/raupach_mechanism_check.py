"""Is the Raupach roughness lever actually connected?  A one-year paired check.

WHY THIS EXISTS.  The 50-year arm 11H is scored on Siberian DJF soil and T2m, whose
detection thresholds are 30-year numbers (+-0.37 K soil, +-0.61 K T2m).  Those say
nothing at all at one year.  So a long run is the only way to see the RESPONSE -- and a
terrible way to discover that the lever was never connected in the first place.  A
coupling field can be computed, sent, received and then quietly ignored, and every
symptom of that is indistinguishable from "the physics did not do much".

This scores the MECHANISM instead.  Roughness enters the momentum budget directly, so
surface stress responds in the first timestep, not the thirtieth year.

THE PAIR.  Tuning_test_11H0_raupach_null (OFF) and Tuning_test_11H1_raupach_1yr (ON),
same binary, same 8-entry namcouple, same initial state.  Only ifraupachz0 and
ECE_CPL_LPJG_Z0 differ.

CORRECTED 2026-08-12.  This used to say the pair required the UNGATED build (11a7debf),
"because with the switch off the gated one stops putting the field while OpenIFS still
expects it".  That was wrong, and reading the source settles it: def_var is
UNCONDITIONAL at OasisCoupler.cpp:303 and only the put at :529 is gated.  So the gated
build declares GUE_Z0HV either way -- which is what OASIS's "namcouple variable not
used" check at mod_oasis_coupler.F90:1266 actually tests -- and OpenIFS's get side is
gated by ECE_CPL_LPJG_Z0 in the same run.  The gated build (f6873559) therefore forms
the pair correctly AND is the build that would go to production, so it is the right one
to test.  Both arms must be on it; 11H0's first year was run on 11a7debf and should be
repeated on f6873559 before the pair is read as strictly one-variable.

THREE DEFECTS had to be cleared before the ON arm executed a single timestep, and each
hid the next:
  1. ECE_CPL_LPJG_Z0 was declared in surfece.F90's NAMECECFG but not in ecearth.F90's.
     NAMECECFG is read from the same fort.4 by both, and a Fortran namelist read aborts
     on any name the reading module does not declare -- forrtl severe (19) on every rank.
  2. bin/guess held the production build 8c5ab467, which contains no Raupach code at
     all, so GUE_Z0HV was never def_var'd and OASIS aborted on the unused namcouple
     entry.  The comp- script builds into lpj_guess*/build/ and does not install, so
     bin/ keeps whatever was last copied there.
  3. couple_put called put_2d(GUE_Z0HV) unconditionally, which aborts whenever the
     namcouple lacks the field -- this is what killed 11G's legs 3-5, not the pair.

WHAT IS PREDICTED, from the offline calculation in plot_raupach_z0.py:
    DJF grid-box z0 over Siberia   0.058 -> 0.140 m   (2.4x)   stress UP
    JJA grid-box z0                0.190 -> 0.120 m   (0.6x)   stress DOWN, smaller
So the DJF and JJA signals should have OPPOSITE signs.  That is a much stronger test
than either one alone -- a plumbing error or a stuck field gives zero in both seasons,
and a sign error gives the wrong sign in both.

THRESHOLDS FIRST, as the campaign requires.  Each metric's 1-year detection threshold is
computed from 11G's own interannual scatter over its completed years, as
1.96*sd*sqrt(2) for a pair of single years.  A difference inside that band is not a
result, however suggestive it looks.

A TRAP THIS AVOIDS.  ewss/nsss are ACCUMULATED turbulent stresses (N m-2 s over the
output step), exactly like the radiation fields, so they must be divided by the
accumulation period.  Forgetting that is the error that made an early M-series pass come
out 3600x too large.
"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

RT = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
OFF = f'{RT}/Tuning_test_11H0_raupach_null'
ON = f'{RT}/Tuning_test_11H1_raupach_1yr'
CTL = f'{RT}/Tuning_test_11G_inppmin50k'          # for the interannual scatter
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/'
        'Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/outdata/oifs/'
        'atm_remapped_1m_lsm_1350-1350.nc')
YEAR = 1350
ACC = 3600.0                       # accumulation period for ewss/nsss
SIB = (55.0, 75.0, 60.0, 180.0)
SEASONS = {'DJF': [11, 0, 1], 'MAM': [2, 3, 4], 'JJA': [5, 6, 7], 'SON': [8, 9, 10]}

print(__doc__)
print('=' * 100)

with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values)
    lat, lon = d['lat'].values, d['lon'].values
if lsm.ndim == 3:
    lsm = lsm[0]
box = np.zeros_like(lsm, bool)
box[np.ix_((lat >= SIB[0]) & (lat < SIB[1]), (lon >= SIB[2]) & (lon <= SIB[3]))] = True
box &= lsm > 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], lsm.shape)


def load(root, var, year):
    """Monthly field for one year, or None.  Accumulated fields are de-accumulated."""
    f = f'{root}/outdata/oifs/atm_remapped_1m_{var}_{year}-{year}.nc'
    if not os.path.exists(f):
        return None
    with xr.open_dataset(f, decode_times=False) as d:
        a = np.squeeze(d[var].values)
    if a.shape[0] != 12:
        return None
    return a / ACC if var in ('ewss', 'nsss', 'sshf', 'slhf') else a


def metrics(root, year):
    """The four mechanism diagnostics, as monthly Siberian box means."""
    ew, ns = load(root, 'ewss', year), load(root, 'nsss', year)
    u, v = load(root, '10u', year), load(root, '10v', year)
    t2, sk = load(root, '2t', year), load(root, 'skt', year)
    if any(x is None for x in (ew, ns, u, v, t2, sk)):
        return None
    out = {}
    fields = {'|tau| [N/m2]': np.sqrt(ew ** 2 + ns ** 2),
              '10m wind [m/s]': np.sqrt(u ** 2 + v ** 2),
              '2t-skt [K]': t2 - sk,
              '2t [K]': t2}
    for name, fld in fields.items():
        out[name] = np.array([float(np.average(fld[k][box], weights=W[box]))
                              for k in range(12)])
    return out


# ---------------------------------------------------------------- thresholds first
print('\n1-YEAR DETECTION THRESHOLDS, from 11G interannual scatter (1.96*sd*sqrt(2))\n')
years = [y for y in range(1350, 1380)
         if os.path.exists(f'{CTL}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc')]
scatter = [metrics(CTL, y) for y in years]
scatter = [s for s in scatter if s]
thr = {}
if len(scatter) >= 5:
    print(f'  {"metric":16s}' + ''.join(f'{s:>12s}' for s in SEASONS))
    for name in scatter[0]:
        thr[name] = {}
        row = f'  {name:16s}'
        for s, mo in SEASONS.items():
            vals = np.array([np.mean(sc[name][mo]) for sc in scatter])
            thr[name][s] = 1.96 * vals.std(ddof=1) * np.sqrt(2)
            row += f'{thr[name][s]:12.4f}'
        print(row)
    print(f'\n  (from {len(scatter)} control years: {years[0]}-{years[-1]})')
else:
    print(f'  only {len(scatter)} control years available -- thresholds not computable')

# ---------------------------------------------------------------- the pair
off, on = metrics(OFF, YEAR), metrics(ON, YEAR)
if off is None or on is None:
    print(f'\n  OFF arm present: {off is not None}   ON arm present: {on is not None}')
    print('  Both arms must have completed year %d.  Nothing to compare yet.' % YEAR)
    sys.exit(0)

print('\n' + '=' * 100)
print(f'\nON minus OFF, Siberian box, year {YEAR}.  * = outside the 1-year threshold\n')
print(f'  {"metric":16s}' + ''.join(f'{s:>14s}' for s in SEASONS))
verdict = {}
for name in off:
    row = f'  {name:16s}'
    for s, mo in SEASONS.items():
        d = np.mean(on[name][mo]) - np.mean(off[name][mo])
        t = thr.get(name, {}).get(s)
        mark = '*' if (t and abs(d) > t) else ' '
        verdict[(name, s)] = (d, t, bool(t and abs(d) > t))
        row += f'{d:+13.4f}{mark}'
    print(row)

# ---------------------------------------------------------------- the verdict
print('\n' + '=' * 100)
tau_djf, _, tau_djf_sig = verdict[('|tau| [N/m2]', 'DJF')]
tau_jja, _, tau_jja_sig = verdict[('|tau| [N/m2]', 'JJA')]
print('\nPRE-REGISTERED READING\n')
print(f'  DJF |tau| {tau_djf:+.4f} ({"resolved" if tau_djf_sig else "NOT resolved"}), '
      f'predicted UP   (z0 0.058 -> 0.140 m)')
print(f'  JJA |tau| {tau_jja:+.4f} ({"resolved" if tau_jja_sig else "NOT resolved"}), '
      f'predicted DOWN (z0 0.190 -> 0.120 m)')
if tau_djf_sig and tau_djf > 0 and tau_jja < 0:
    print('\n  -> LEVER CONNECTED, and the seasonal sign reversal matches the prediction.')
    print('     Spend the 50 years: 11H is worth running.')
elif not tau_djf_sig:
    print('\n  -> NOT CONNECTED, or too weak to see.  The field is being sent and')
    print('     apparently ignored.  Check the IFS side (SURFECE_GET_Z0VEG, vupdz0)')
    print('     before spending anything on a long run.')
else:
    print('\n  -> MOVES, BUT NOT AS PREDICTED.  A sign error in vupdz0 or in the')
    print('     drag-space aggregation is the first thing to check.')
print()
