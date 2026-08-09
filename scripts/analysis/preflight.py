"""preflight.py PARAM [PARAM ...] -- answer "has this been tried, and what is it?"

WHY THIS EXISTS.  Of the 38 claims this campaign has had to retract, roughly 45 % were
answerable at a desk in seconds and were not checked.  Six 48-year runs were spent on
them.  The three failure modes this tool removes, with the run each one cost:

  * PROPOSING SOMETHING ALREADY RUN.  RCL_OVERLAPLIQICE 0.65->0.35 was submitted as a
    new lever on 2026-08-09.  It is A1b, run for 44 years in round 3.  The A-series
    label "ovl" means the liquid/ice DEPOSITION overlap; it was read as the radiative
    cloud overlap.  A label is not an identifier -- the namelist value is.

  * READING A SOURCE DEFAULT AS THE RUNTIME VALUE.  On the same day, "NAERCLD=0, so the
    aerosol-cloud path is off and CCN is a uniform 125" came from a default assignment
    in sucldp.F90.  The runtime fort.4 carries LMACV2SP_CCNF=.true., i.e. a second CCN
    path is active, and PCCN is a field.  Defaults are what applies when the namelist is
    SILENT; only fort.4 says what the model ran with.

  * NOT ASKING WHY A VALUE IS WHAT IT IS.  Lowering RVICE was proposed as a
    tropics-warming lever without checking that RVICE had been RAISED from its 0.13
    source default to 0.16 specifically to suppress Southern Ocean high cloud.  That is
    in the project issues, not in the code.

WHAT IT REPORTS, per parameter:
  1. the source default, and the file/line that sets it
  2. whether it is reachable from a namelist at all (and which one)
  3. every run that has ever SET it, with the value, straight from each run's fort.4
  4. every runscript that mentions it, including withdrawn ones
  5. git log touching it, so a non-default value has provenance
  6. every mention in the report and notes, so prior findings surface

It reads the model's own files; it does not consult any cached list that could go stale.

USAGE
    python3 preflight.py RVICE
    python3 preflight.py RCL_OVERLAPLIQICE RCLDIFF RCCN
"""
import os, re, subprocess, sys, glob

BASE = ('/work/ab0246/a270092/postprocessing/'
        'investigation_awiesm3_high_lat_cold_bias_global_toa_positive')
SRC = '/work/ab0246/a270092/model_codes/oifsamip-cy48/oifs-48r1/ifs-source'
RUNTIMES = ['/work/bb1469/a270092/runtime/oifsamip-cy48',
            '/work/bb1469/a270270/runtime/awiesm3-v3.4',
            '/work/bb1469/a270092/runtime/awiesm3-v3.4']
SCRIPTS = '/home/a/a270092/esm_tools/runscripts/oifsamip'


def sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                              timeout=120).stdout
    except Exception:
        return ''


def default_of(p):
    """The source default and where it is set."""
    out = sh(f"grep -rn --include=*.F90 -E '^ *{p} *=[^=]' {SRC} | head -5")
    return [l for l in out.splitlines() if l.strip()]


def namelists_with(p):
    hits = []
    for f in glob.glob(f'{SRC}/**/*.nam.h', recursive=True):
        try:
            txt = open(f).read()
        except OSError:
            continue
        if re.search(rf'\b{p}\b', txt):
            m = re.search(r'NAMELIST */([A-Z0-9_]+)/', txt)
            hits.append(m.group(1) if m else os.path.basename(f))
    return hits


def runs_setting(p):
    """Every run whose fort.4 sets this parameter -- the authoritative record."""
    found = []
    for rt in RUNTIMES:
        if not os.path.isdir(rt):
            continue
        for d in sorted(os.listdir(rt)):
            fs = sorted(glob.glob(f'{rt}/{d}/run_*/work/fort.4'), reverse=True)
            if not fs:
                continue
            try:
                txt = open(fs[0]).read()
            except OSError:
                continue
            m = re.search(rf'^\s*{p}\s*=\s*(\S+)', txt, re.M | re.I)
            if m:
                found.append((d, m.group(1).rstrip(',')))
    return found


def main(params):
    for p in params:
        P = p.upper()
        print('=' * 94)
        print(f'PREFLIGHT: {P}')
        print('=' * 94)

        print('\n1. SOURCE DEFAULT (applies only when the namelist is silent)')
        d = default_of(P)
        print('\n'.join('   ' + x for x in d) if d else '   no plain assignment found')

        print('\n2. NAMELIST REACHABILITY')
        nl = namelists_with(P)
        if nl:
            print(f'   in {", ".join(sorted(set(nl)))} -- settable without a rebuild')
        else:
            print('   NOT in any namelist. Setting it needs a source change '
                  '(declaration + association + namelist entry).')

        print('\n3. RUNS THAT SET IT  <- the authoritative "has this been tried"')
        rs = runs_setting(P)
        if rs:
            vals = {}
            for name, v in rs:
                vals.setdefault(v, []).append(name)
            for v, names in sorted(vals.items()):
                short = [n.replace('Tuning_test_', '').replace('amip_', '') for n in names]
                print(f'   = {v:<12s} {len(names):3d} run(s): '
                      f'{", ".join(short[:8])}{" ..." if len(short) > 8 else ""}')
            print(f'   *** ALREADY EXERCISED at {len(vals)} distinct value(s). '
                  f'Check these before proposing a new one.')
        else:
            print('   never set in any run -- genuinely untested')

        print('\n4. RUNSCRIPTS MENTIONING IT (including withdrawn)')
        o = sh(f"grep -rl '{P}' {SCRIPTS} 2>/dev/null | head -12")
        for l in o.splitlines():
            tag = ' [WITHDRAWN]' if '/withdrawn/' in l else ''
            print(f'   {os.path.basename(l)}{tag}')
        if not o.strip():
            print('   none')

        print('\n5. PROVENANCE (git log in the model tree)')
        o = sh(f"cd {SRC} && git log --oneline -S'{P}' -- . 2>/dev/null | head -6")
        print('\n'.join('   ' + x for x in o.splitlines()) if o.strip()
              else '   no commits touch it (or not a git tree)')

        print('\n6. PRIOR FINDINGS in the report and notes')
        o = sh(f"grep -rn '{P}' {BASE}/report/report.tex {BASE}/notes/*.md 2>/dev/null "
               f"| head -8")
        for l in o.splitlines():
            f, n, rest = (l.split(':', 2) + ['', ''])[:3]
            print(f'   {os.path.basename(f)}:{n}  {rest.strip()[:88]}')
        if not o.strip():
            print('   not mentioned -- no prior finding to contradict')
        print()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
