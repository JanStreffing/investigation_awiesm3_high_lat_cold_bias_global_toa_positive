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

  * LOOKING ONLY AT THIS CAMPAIGN.  Fixed 2026-08-10, after the tool said "never set in
    any run" TWICE IN ONE DAY about parameters with real AWI history.  RPRCON had been
    tuned to 0.7E-3 at TCO319 (project_management #87) and launched once at TCO95 (#95);
    RSNOWLIN2 had been set to 0.04 in the awiesm3 v3.4.2 TCO95L91-DARS2 runscripts and to
    0.025 in awicm3 MELDPOND -- along with a whole jointly-tuned cloud stack (DETRPEN
    1.32E-4, ENTRORG 2.07E-3, RMFDEPS 0.48, ENTRDD 1.08E-4, RVICE 0.18, RLCRITSNOW
    1.46E-5).  The tool was scanning ONE runtime tree and ONE runscript directory.  It now
    scans every runtime tree and every runscript directory on disk, and reports the VALUE,
    because "mentioned in a runscript" and "set to something in a runscript" are different
    facts and only the second one is prior art.

  WHAT IT STILL CANNOT SEE: the GitHub issues.  Sections 3-5 cover disk only.  For any
  parameter that matters, also search AWI-ESM/project_management and
  tsemmler05/AWI-CM3-HighResMIP -- that is where the WHY lives, and issue search does not
  index comment bodies, so fetch the comments and grep them.

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

# EVERY runtime tree and EVERY runscript directory, not just this campaign's.  These are
# globs on purpose: a hard-coded list is what let RPRCON and RSNOWLIN2 read as "never
# set" on 2026-08-10 when both had AWI history one directory over.
RUNTIME_ROOTS = ['/work/bb1469/a270092/runtime', '/work/bb1469/a270270/runtime',
                 '/work/ab0246/a270092/runtime', '/work/bb1469/a270092/runtime_old']
SCRIPT_ROOTS = ['/home/a/a270092/esm_tools/runscripts']
# The campaign's own setup, so output can separate "tried here" from "tried at AWI".
HOME_SETUP = 'oifsamip'


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
    """Every run on disk whose fort.4 sets this, across ALL setups.

    Returns (setup, experiment, value).  fort.4 is the authority on what a run actually
    used -- runscripts can be edited after the fact, fort.4 cannot.
    """
    found = []
    for root in RUNTIME_ROOTS:
        if not os.path.isdir(root):
            continue
        for setup in sorted(os.listdir(root)):
            sd = f'{root}/{setup}'
            if not os.path.isdir(sd):
                continue
            for d in sorted(os.listdir(sd)):
                fs = sorted(glob.glob(f'{sd}/{d}/run_*/work/fort.4'), reverse=True)
                if not fs:
                    continue
                try:
                    txt = open(fs[0]).read()
                except OSError:
                    continue
                m = re.search(rf'^\s*{p}\s*=\s*(\S+)', txt, re.M | re.I)
                if m:
                    found.append((setup, d, m.group(1).rstrip(',')))
    return found


def scripts_setting(p):
    """Every runscript on disk that ASSIGNS this parameter, with its value and setup.

    Distinguishes assignment (`RSNOWLIN2: 0.04`) from mere mention (a comment), because
    only the first is prior art.  This is the check that was missing.
    """
    found = []
    for root in SCRIPT_ROOTS:
        for f in glob.glob(f'{root}/**/*.y*ml', recursive=True):
            try:
                txt = open(f).read()
            except OSError:
                continue
            m = re.search(rf'^\s*["\']?{p}["\']?\s*:\s*([^\s#]+)', txt, re.M | re.I)
            rel = os.path.relpath(f, root)
            setup = rel.split(os.sep)[0]
            if m:
                found.append((setup, rel, m.group(1).rstrip(',')))
            elif re.search(rf'\b{p}\b', txt, re.I):
                found.append((setup, rel, None))          # mentioned only
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

        print('\n3. RUNS THAT SET IT, ALL SETUPS  <- what the model actually ran with')
        rs = runs_setting(P)
        if rs:
            vals = {}
            for setup, name, v in rs:
                vals.setdefault(v, []).append((setup, name))
            for v, items in sorted(vals.items()):
                here = sum(1 for s, _ in items if s.startswith(HOME_SETUP))
                short = [n.replace('Tuning_test_', '').replace('amip_', '')
                         for _, n in items]
                where = f'{here} here, {len(items)-here} elsewhere' if here != len(items) \
                        else 'this campaign'
                print(f'   = {v:<12s} {len(items):3d} run(s) [{where}]: '
                      f'{", ".join(short[:6])}{" ..." if len(short) > 6 else ""}')
            print(f'   *** ALREADY EXERCISED at {len(vals)} distinct value(s). '
                  f'Check these before proposing a new one.')
        else:
            print('   no fort.4 on disk sets it')

        print('\n4. RUNSCRIPTS THAT SET IT, ALL SETUPS  <- prior art even if never run')
        ss = scripts_setting(P)
        assigned = [x for x in ss if x[2] is not None]
        mentioned = [x for x in ss if x[2] is None]
        if assigned:
            for setup, rel, v in sorted(assigned, key=lambda t: (t[0], t[1])):
                tag = ' [WITHDRAWN]' if '/withdrawn/' in rel else ''
                mark = '' if setup.startswith(HOME_SETUP) else '  <- OTHER SETUP'
                print(f'   = {v:<12s} {setup}/{os.path.basename(rel)}{tag}{mark}')
            if any(not s.startswith(HOME_SETUP) for s, _, _ in assigned):
                print('   *** SET OUTSIDE THIS CAMPAIGN. That is prior art: find out what '
                      'it was for\n       before proposing a value.')
        if mentioned:
            print(f'   mentioned only (no assignment) in {len(mentioned)} file(s): '
                  f'{", ".join(os.path.basename(r) for _, r, _ in mentioned[:5])}')
        if not ss:
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
