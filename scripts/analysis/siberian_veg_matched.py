"""Siberian vegetation, compared on a MATCHED decade rather than "last on disk".

WHY THIS EXISTS.  coupled_11IJ_deep.py section 4 reads the last decade each arm has,
which is correct when every arm is the same length and wrong the moment they are not.
Scoring 11L/11M (40 yr, last decade 1380-89) against 11G (50 yr, last decade 1390-99)
compares different decades of a DECLINING forest, so part of any difference is just
drift.  On 2026-08-19 that inflated the apparent damage; this script removes it.

WHAT IT REPORTS.  Area-weighted Siberian box mean of the final year each cell reports in
the chosen leg, from fpc.out -- written per rank by LPJ-GUESS as the run proceeds, so it
survives a finalisation abort.  Handles both output layouts: outdata/lpj_guess/<leg>/run*/
and the in-work run*/output/ tree that short runs leave behind.

BNS is boreal needleleaf summergreen (larch), the east-Siberian dominant.  Standing
finding, re-derived repeatedly: BNS is NOT climate-gate limited, it is outcompeted by C3
grass.  So read a BNS change as competition or productivity, never as a gate opening.

Usage:  siberian_veg_matched.py <leg> <tag>=<dir> [<tag>=<dir> ...]
        the FIRST arm given is the baseline everything else is differenced against.
"""
import glob
import sys
import numpy as np
import warnings
warnings.filterwarnings('ignore')

SIB = (55.0, 75.0, 60.0, 180.0)
COLS = ['BNE', 'BINE', 'BNS', 'IBS', 'TREEFPC', 'C3G', 'GRASSFPC', 'Total']


def fpc(root, leg):
    files = sorted(glob.glob(f'{root}/outdata/lpj_guess/{leg}/run*/fpc.out'))
    if not files:
        files = sorted(glob.glob(f'{root}/run_{leg}/work/run*/output/fpc.out'))
    if not files:
        return None, 0
    tot = {c: 0.0 for c in COLS}
    W = 0.0
    n = 0
    for fn in files:
        try:
            with open(fn) as fh:
                hdr = fh.readline().split()
                if 'Lat' not in hdr:
                    continue
                idx = {c: hdr.index(c) for c in COLS if c in hdr}
                ila, ilo, iyr = hdr.index('Lat'), hdr.index('Lon'), hdr.index('Year')
                last = {}
                for ln in fh:
                    p = ln.split()
                    if len(p) < len(hdr):
                        continue
                    key = (p[ilo], p[ila])
                    if key not in last or int(p[iyr]) > int(last[key][iyr]):
                        last[key] = p
                for p in last.values():
                    la, lo = float(p[ila]), float(p[ilo]) % 360
                    if not (SIB[0] <= la <= SIB[1] and SIB[2] <= lo <= SIB[3]):
                        continue
                    w = np.cos(np.deg2rad(la))
                    W += w
                    n += 1
                    for c in COLS:
                        tot[c] += float(p[idx[c]]) * w if c in idx else 0.0
        except OSError:
            pass
    if W == 0:
        return None, 0
    return {c: tot[c] / W for c in COLS}, n


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    leg = sys.argv[1]
    arms = [a.split('=', 1) for a in sys.argv[2:]]
    print(__doc__.split('Usage:')[0])
    print(f'Siberian box, MATCHED decade {leg} (final year each cell reports)\n')
    print(f'{"arm":6s} {"cells":>6s}  ' + '  '.join(f'{c:>8s}' for c in COLS))
    base = None
    for tag, root in arms:
        v, n = fpc(root, leg)
        if v is None:
            print(f'{tag:6s} no fpc.out for this leg')
            continue
        print(f'{tag:6s} {n:6d}  ' + '  '.join(f'{v[c]:8.4f}' for c in COLS))
        if base is None:
            base = v
        else:
            print(f'{"":6s} {"vs base":>6s}  '
                  + '  '.join(f'{v[c] - base[c]:+8.4f}' for c in COLS))
            print(f'{"":6s} {"%":>6s}  ' + '  '.join(
                f'{100 * (v[c] - base[c]) / base[c]:+7.1f}%' if base[c] else '      --'
                for c in COLS))


if __name__ == '__main__':
    main()
