"""Score the restart_target_continuity A/B, with the instrument checked before the data.

THE TEST.  Two runs, one binary, one namelist line apart, both fixed_LU:

    <ctl arm>   restart_target_continuity 0    pre-fix: frac_old re-summed from stands
    <fix arm>   restart_target_continuity 1    frac_old differenced against the carried value

Both carry print_lc_change_diag 1, which makes lc_changed() emit one LCDIAG line per
gridcell per year whenever the land-cover machinery sees any nonzero change.

PRE-REGISTERED CRITERION (notes/lpjg_merge_plan.md, condition 2.1):

    ctl LCDIAG > 0  AND  fix LCDIAG == 0   ->  PASS
    ctl LCDIAG == 0                        ->  NO TEST, redesign; explicitly NOT a pass

WHY THE CONTROL IS CHECKED FIRST.  Four instruments in a row returned a plausible zero
on this measurement before any of them was actually measuring:

  1. Slurm job state.  esm_tools reports COMPLETED on a crashed leg and advances the
     date file, so a 72-second death scored as a finished run.  -> score on output.
  2. The log path.  LPJ-GUESS dprintf lands in run_*/work/**/guess*.log, NOT log/*.log.
  3. Variable ordering.  The LCDIAG call sat ~40 lines above where change_st is
     accumulated and before change_gross_lcc existed, so it reported change_st as 0.0
     unconditionally -- the stand-type term, the granularity the defect lives at.
     Fixed in lpjg 4d91c3b; needs binary guess.restart_ab_lcdiagfix (md5 6022479c).
  4. A typo'd needle.  restart_stress_compare.py grepped 'Carried land-cover metadata
     disagrees'; the string is 'Carried stand-type metadata disagrees'.  Guard-fallback
     counts were 0 by construction.

So a zero from the control arm is treated as an instrument fault until proven otherwise,
and this script will not report PASS unless the control fired.

Reads only per-rank LPJ-GUESS text output, written as the run proceeds, which survives a
finalisation abort.

Usage:
    python lcdiag_ab_score.py <fix_arm_dir> <ctl_arm_dir>
"""
import glob, os, re, sys

# Verified against lpj_guess_repairfirst 2026-08-18.  Do not edit without re-grepping the
# source -- see instrument failure 4 above.
N_LCDIAG   = 'LCDIAG '
N_GUARD    = 'Carried stand-type metadata disagrees'
N_MISMATCH = 'LUH3 previous LC/ST mismatch'
N_FATAL    = ('does not sum to one', 'too large to be input roundoff',
              'exceeds available LUH3', 'failed to read index for state file',
              'Physical restart stand area is not usable')

RE_LCDIAG = re.compile(
    r'LCDIAG year=(\d+) changeLC=(\S+) change_st=(\S+) (?:gross=(\S+) )?fsc=(\d)')
RE_LCDRIFT = re.compile(r'LCDRIFT year=(\d+) maxdrift=(\S+) st=(\d+) limit=(\S+)')


def legs(arm):
    return sorted(glob.glob(os.path.join(arm, 'run_*')))


def scan(leg):
    o = dict(lcdiag=0, guard=0, mismatch=0, fatal=[], legacy=False,
             max_lc=0.0, max_st=0.0, fsc=0, drift=[], limit=None)
    for fn in glob.glob(os.path.join(leg, 'work', '**', 'guess*.log'), recursive=True):
        try:
            with open(fn, errors='ignore') as fh:
                for ln in fh:
                    if N_LCDIAG in ln:
                        o['lcdiag'] += 1
                        m = RE_LCDIAG.search(ln)
                        if m:
                            if m.group(4) is None:
                                o['legacy'] = True   # binary predates 4d91c3b
                            o['max_lc'] = max(o['max_lc'], abs(float(m.group(2))))
                            o['max_st'] = max(o['max_st'], abs(float(m.group(3))))
                            o['fsc'] += int(m.group(5))
                    if 'LCDRIFT ' in ln:
                        m = RE_LCDRIFT.search(ln)
                        if m:
                            o['drift'].append(float(m.group(2)))
                            o['limit'] = float(m.group(4))
                    if N_GUARD in ln:
                        o['guard'] += 1
                    if N_MISMATCH in ln:
                        o['mismatch'] += 1
                    if any(p in ln for p in N_FATAL):
                        o['fatal'].append(ln.strip()[:110])
        except OSError:
            pass
    return o


def years(leg):
    """Non-empty fpc.out.  A leg that never finished a year writes them empty; this is
    the progress signal, NOT the Slurm job state."""
    return sum(1 for f in glob.glob(os.path.join(leg, 'work', '**', 'fpc.out'),
                                    recursive=True) if os.path.getsize(f) > 0)


def report(arm):
    L = legs(arm)
    print('\n%s   %d leg(s)' % (os.path.basename(arm.rstrip('/')), len(L)))
    if not L:
        print('   no legs yet')
        return None
    t = dict(lcdiag=0, guard=0, mismatch=0, fatal=0, years=0, legacy=False,
             nlegs=len(L), guard_after_leg1=0, limit=None)
    alldrift = []
    print('   %24s %4s %8s %8s %6s %6s %11s'
          % ('leg', 'yrs', 'LCDIAG', 'LCDRIFT', 'guard', 'LC/ST', 'maxdrift'))
    for i, leg in enumerate(L):
        s = scan(leg)
        y = years(leg)
        for k in ('lcdiag', 'guard', 'mismatch'):
            t[k] += s[k]
        t['fatal'] += len(s['fatal'])
        t['years'] += y
        t['legacy'] = t['legacy'] or s['legacy']
        if i > 0:
            t['guard_after_leg1'] += s['guard']
        md = max(s['drift']) if s['drift'] else 0.0
        alldrift.extend(s['drift'])
        if s['limit'] is not None:
            t['limit'] = s['limit']
        print('   %24s %4d %8d %8d %6d %6d %11.3e'
              % (os.path.basename(leg), y, s['lcdiag'], len(s['drift']), s['guard'],
                 s['mismatch'], md))
        for f in s['fatal'][:2]:
            print('        FATAL: %s' % f)
    if t['legacy']:
        print('   !! binary predates lpjg 4d91c3b: change_st was read before it was')
        print('      computed, so these counts are NOT scorable.')
    # Is LC_PREV_FRAC_MAX_DRIFT anywhere near the data it guards?  Zero fallbacks only
    # means something if the bound is comparable to the observed spread; if the drift is
    # ~1e-7 against a 0.1 limit, the guard is decorative and "no fallbacks" is not
    # evidence that the carried value was validated.
    if alldrift:
        alldrift.sort()
        n = len(alldrift)
        lim = t['limit']
        print('   carried-vs-physical drift: n=%d  median=%.3e  p99=%.3e  max=%.3e'
              % (n, alldrift[n // 2], alldrift[min(n - 1, int(0.99 * n))], alldrift[-1]))
        if lim:
            print('   limit=%.3g  -> max observed is %.3g x the limit'
                  % (lim, alldrift[-1] / lim))
            if alldrift[-1] < lim / 1000.0:
                print('   NOTE: every observed drift is >1000x below the limit. Zero')
                print('         fallbacks is therefore not evidence the guard works --')
                print('         it has never been near firing. Condition 5 is untested.')
    t['drift_max'] = alldrift[-1] if alldrift else None
    return t


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    print(__doc__.split('Usage:')[0])
    print('=' * 96)
    fix = report(sys.argv[1])
    ctl = report(sys.argv[2])
    print('\n' + '=' * 96)
    if not fix or not ctl:
        print('VERDICT: incomplete -- one arm has no legs.'); return
    if fix['legacy'] or ctl['legacy']:
        print('VERDICT: NOT SCORABLE. Re-run on guess.restart_ab_lcdiagfix (6022479c).')
        return
    if fix['years'] == 0 or ctl['years'] == 0:
        print('VERDICT: NO TEST. An arm never completed a model year, so LCDIAG was')
        print('         never reached. Zero here means nothing.'); return
    if fix['nlegs'] < 2 or ctl['nlegs'] < 2:
        print('VERDICT: NO TEST. Needs >=2 legs; the restart boundary is the point.')
        return
    if fix['fatal'] or ctl['fatal']:
        print('VERDICT: FAIL. Hard failures present; fix them before scoring.'); return

    print('   control arm LCDIAG total : %d' % ctl['lcdiag'])
    print('   fix arm     LCDIAG total : %d' % fix['lcdiag'])
    # The guard sits behind `if(use_prev)` and use_prev = restart_target_continuity, so
    # it NEVER runs in the control arm.  A ctl guard count of 0 is uninformative by
    # construction; merge condition 5 is testable in the fix arm only.
    print('   guard fallbacks after leg 1 : fix %d   (expected 0; leg 1 is a foreign '
          'parent)' % fix['guard_after_leg1'])
    print('                                 ctl n/a -- guard is gated off when '
          'restart_target_continuity 0')
    print('   LUH3 previous LC/ST mismatch : fix %d, ctl %d   (expected 0)'
          % (fix['mismatch'], ctl['mismatch']))
    print()

    if ctl['lcdiag'] == 0:
        print('VERDICT: NO TEST -- and suspect the instrument first.')
        print('   The control reproduces the pre-fix behaviour by construction, so a zero')
        print('   from it means the diagnostic did not fire, not that nothing moved.')
        print('   Four instruments have already failed this way. NOT a pass.')
        return
    if fix['lcdiag'] > 0:
        print('VERDICT: FAIL. The fix arm still registers %d changes under fixed_LU,'
              % fix['lcdiag'])
        print('   where differencing target against target must give exactly zero.')
        return
    print('VERDICT: PASS. Control fired (%d), fix is exactly zero across every leg.'
          % ctl['lcdiag'])
    if fix['guard_after_leg1']:
        print('   caveat: guard fallbacks after leg 1 -- LC_PREV_FRAC_MAX_DRIFT mis-tuned')
    if fix['mismatch'] or ctl['mismatch']:
        print('   caveat: LUH3 previous LC/ST mismatch present')


if __name__ == '__main__':
    main()
