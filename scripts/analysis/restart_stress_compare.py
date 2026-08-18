"""Did differencing against the carried value stop the leg-boundary wobble?

THE TEST.  Four 3-leg runs from the same 1350 state, differing in two things only:

    R1_pi_ctl   fixed_LU 1850    guess.peatlu_policies            (pre-fix)
    R2_pi_new   fixed_LU 1850    guess.restart_targetcontinuity   (lpjg ab45fdd)
    R3_tr_ctl   transient LU     guess.peatlu_policies            (pre-fix)
    R4_tr_new   transient LU     guess.restart_targetcontinuity

Three 1-year legs give TWO restart boundaries.  Leg 1 reads the 2000-year offline
spin-up -- a foreign parent configuration, the case getlandcover's original comment was
written for -- so the leg1->leg2 boundary is contaminated by construction and is not the
evidence.  **The leg2->leg3 boundary is the only clean self-restart** and is what this
script scores.

WHAT PASS LOOKS LIKE.
  fixed_LU (R2 vs R1): crp/pst and the peat-cell count IDENTICAL across leg2->leg3.
      The LUH3 target is the same every year, so differencing target-to-target must give
      frac_change == 0 exactly.  R1 is expected to wobble; that is the control.
  transient (R4 vs R3): land use genuinely changes 1350->1352, so crp/pst SHOULD move.
      What must not happen is a crash, a guard fallback storm, or peat blinking.
      R4 completing at all is also the first exercise of peat_lu_conflict_policy 0's
      transfer path -- under fixed_LU that path refuses instead.

Reads only per-leg LPJ-GUESS text output, which is written per rank as the run proceeds
and therefore survives a finalisation abort.
"""
import os, glob, sys
import numpy as np

B = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
RUNS = [('R1_pi_ctl', 'fixed_LU', 'pre-fix'),
        ('R2_pi_new', 'fixed_LU', 'restart fix'),
        ('R3_tr_ctl', 'transient', 'pre-fix'),
        ('R4_tr_new', 'transient', 'restart fix')]
LUT = ('crp', 'pst', 'psl')


def legs(run):
    """Leg directories, newest layout first.

    These short runs leave the annual .out files in the WORK tree under
    run<N>/output/, and only the daily ones are moved to outdata/lpj_guess/<leg>/run<N>/.
    Longer campaign runs (11I) have them directly under outdata. Handle both.
    """
    out = sorted(glob.glob(f'{B}/Restart_{run}/run_*/work'))
    if out:
        return out
    return sorted(glob.glob(f'{B}/Restart_{run}/outdata/lpj_guess/*'))


def _find(leg, name):
    hits = sorted(glob.glob(f'{leg}/run*/output/{name}'))
    if not hits:
        hits = sorted(glob.glob(f'{leg}/run*/{name}'))
    return hits


def lut_means(leg):
    """Area-weighted global mean of the land-use tile fractions, last year in the leg."""
    per = {}
    for fn in _find(leg, 'fracLut_yearly.out'):
        with open(fn) as fh:
            hdr = fh.readline().split()
            if not all(c in hdr for c in LUT):
                return None, None
            ix = {c: hdr.index(c) for c in LUT}
            ilat = hdr.index('Lat')
            for line in fh:
                p = line.split()
                if len(p) < len(hdr):
                    continue
                per[(p[0], p[1])] = (int(p[2]), float(p[ilat]),
                                     [float(p[ix[c]]) for c in LUT])
    if not per:
        return None, None
    yr = max(v[0] for v in per.values())
    keep = [v for v in per.values() if v[0] == yr]
    w = np.cos(np.deg2rad(np.array([v[1] for v in keep])))
    vals = np.array([v[2] for v in keep])
    return yr, [float(np.average(vals[:, i], weights=w)) for i in range(len(LUT))]


def peat(leg):
    """(year, cells with peat > 0, area-weighted mean PEAT_STANDFRAC)."""
    per = {}
    for fn in _find(leg, 'fpc.out'):
        with open(fn) as fh:
            hdr = fh.readline().split()
            if 'PEAT_STANDFRAC' not in hdr:
                return None, None, None
            ip, ilat = hdr.index('PEAT_STANDFRAC'), hdr.index('Lat')
            for line in fh:
                p = line.split()
                if len(p) < len(hdr):
                    continue
                per[(p[0], p[1])] = (int(p[2]), float(p[ilat]), float(p[ip]))
    if not per:
        return None, None, None
    yr = max(v[0] for v in per.values())
    keep = [v for v in per.values() if v[0] == yr]
    w = np.cos(np.deg2rad(np.array([v[1] for v in keep])))
    psf = np.array([v[2] for v in keep])
    return yr, int((psf > 0).sum()), float(np.average(psf, weights=w))


def guard_fallbacks(run):
    n = 0
    for fn in glob.glob(f'{B}/Restart_{run}/run_*/work/run*/guess*.log'):
        try:
            with open(fn, errors='ignore') as fh:
                n += sum(1 for ln in fh if 'Carried stand-type metadata disagrees' in ln)
        except OSError:
            pass
    return n


def failures(run):
    pats = ('exceeds available LUH3', 'not normalized before land-cover',
            'does not match physical stands', 'leaves no room', 'not usable before')
    hits = []
    for fn in glob.glob(f'{B}/Restart_{run}/run_*/work/run*/guess*.log'):
        try:
            with open(fn, errors='ignore') as fh:
                for ln in fh:
                    if any(p in ln for p in pats):
                        hits.append(ln.strip())
                        break
        except OSError:
            pass
    return hits


print(__doc__)
print('=' * 100)

for run, mode, arm in RUNS:
    L = legs(run)
    print(f'\n{run}   [{mode}, {arm}]   {len(L)} leg(s)')
    if not L:
        print('   no output yet')
        continue
    rows = []
    for leg in L:
        yr, lm = lut_means(leg)
        pyr, pn, ppsf = peat(leg)
        tag = leg.split('/run_')[-1].split('/')[0] if '/run_' in leg else os.path.basename(leg)
        rows.append((tag, yr, lm, pn, ppsf))
    print(f'   {"leg":>22s} {"yr":>5s} {"crp":>11s} {"pst":>11s} {"psl":>11s} '
          f'{"peat>0":>7s} {"meanPSF":>9s}')
    for name, yr, lm, pn, ppsf in rows:
        if lm is None:
            print(f'   {name:>22s}  (no fracLut yet)')
            continue
        print(f'   {name:>22s} {yr:5d} {lm[0]:11.6f} {lm[1]:11.6f} {lm[2]:11.6f} '
              f'{pn if pn is not None else -1:7d} {ppsf if ppsf is not None else float("nan"):9.6f}')

    # the clean self-restart boundary
    if len(rows) >= 3 and rows[1][2] and rows[2][2]:
        a, b = rows[1], rows[2]
        d = [b[2][i] - a[2][i] for i in range(len(LUT))]
        dn = (b[3] - a[3]) if (a[3] is not None and b[3] is not None) else None
        exact = all(x == 0.0 for x in d)
        print(f'\n   leg2->leg3 (CLEAN self-restart):')
        print(f'     d(crp) {d[0]:+.9f}   d(pst) {d[1]:+.9f}   d(psl) {d[2]:+.9f}'
              + ('   <-- EXACTLY ZERO' if exact else ''))
        if dn is not None:
            print(f'     d(peat cells) {dn:+d}')
        if mode == 'fixed_LU':
            print('     expected under fixed_LU: exactly zero'
                  + ('  PASS' if exact and dn == 0 else '  FAIL'))
        else:
            print('     transient: movement is expected; watch for crashes/blinking instead')

    print(f'   guard fallbacks (carried value rejected): {guard_fallbacks(run)}')
    f = failures(run)
    print(f'   hard failures: {len(f)}' + (f'  e.g. {f[0][:90]}' if f else ''))

print('\n' + '=' * 100)
print('Leg 1 is a restart from a FOREIGN parent (the 2000-year spin-up); guard fallbacks')
print('there are expected and correct. Fallbacks at leg2->leg3 would not be.')
