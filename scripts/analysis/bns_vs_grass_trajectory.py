"""Larch against grass in the Siberian box, year by year: is BNS holding, or still being taken?

The standing finding is that BNS (boreal needleleaf summergreen, larch) is NOT climate-gate
limited in this model -- every gate is open across the box -- it is outcompeted by C3 grass.
So this reports the competition, not thresholds: who holds the cover, who makes the carbon,
and on how many cells any larch is left at all.

Per model year over the Siberian box (55-75N, 60-180E), area weighted by cos(lat), from the
LPJ-GUESS per-rank .out text files:
  fpc     BNS, C3G, TREEFPC, GRASSFPC, and the share of cells with GRASSFPC > TREEFPC
  dens    BNS stems per m2 -- near zero means establishment never happens,
          non-zero with low fpc means larch establishes and then loses
  cmass   standing carbon, BNS against C3G
  anpp    productivity, BNS against C3G -- the rate at which each is winning
Plus the count of cells still carrying BNS above 0.01 and above 0.05 fpc, which is the
"is any left" question in its bluntest form.

Usage:  bns_vs_grass_trajectory.py <tag>=<run dir>[:<leg>] [<tag>=<run dir>[:<leg>] ...]
        leg defaults to every leg found under outdata/lpj_guess.
"""
import glob, os, sys
import numpy as np, warnings
warnings.filterwarnings('ignore')

SIB = (55.0, 75.0, 60.0, 180.0)


def read(root, leg, fn, cols):
    """{year: {col: (weighted sum, weight, [values])}} over the Siberian box"""
    files = sorted(glob.glob(f'{root}/outdata/lpj_guess/{leg}/run*/{fn}'))
    if not files:
        files = sorted(glob.glob(f'{root}/run_{leg}/work/run*/output/{fn}'))
    out = {}
    for path in files:
        with open(path) as fh:
            hdr = fh.readline().split()
            if 'Lat' not in hdr:
                continue
            ila, ilo, iyr = hdr.index('Lat'), hdr.index('Lon'), hdr.index('Year')
            idx = {c: hdr.index(c) for c in cols if c in hdr}
            for ln in fh:
                p = ln.split()
                if len(p) < len(hdr):
                    continue
                la, lo = float(p[ila]), float(p[ilo]) % 360
                if not (SIB[0] <= la <= SIB[1] and SIB[2] <= lo <= SIB[3]):
                    continue
                y = int(p[iyr]); w = np.cos(np.deg2rad(la))
                d = out.setdefault(y, {'w': 0.0, 'n': 0})
                d['w'] += w; d['n'] += 1
                for c, i in idx.items():
                    v = float(p[i])
                    d[c] = d.get(c, 0.0) + v * w
                    d.setdefault(c + '#', []).append(v)
    return out


def legs(root):
    return sorted(os.path.basename(p) for p in glob.glob(f'{root}/outdata/lpj_guess/*') if os.path.isdir(p))


print(__doc__.split('Usage:')[0])
for spec in sys.argv[1:]:
    tag, path = spec.split('=', 1)
    root, _, leg = path.partition(':')
    for lg in ([leg] if leg else legs(root)):
        F = read(root, lg, 'fpc.out', ['BNS', 'C3G', 'TREEFPC', 'GRASSFPC'])
        D = read(root, lg, 'dens.out', ['BNS'])
        C = read(root, lg, 'cmass.out', ['BNS', 'C3G'])
        A = read(root, lg, 'anpp.out', ['BNS', 'C3G'])
        if not F:
            print(f'{tag} {lg}: no fpc.out'); continue
        print(f'== {tag}  leg {lg}  ({F[max(F)]["n"]} Siberian cells)')
        print(f'  {"year":>6}{"BNS fpc":>9}{"C3G fpc":>9}{"TREE":>8}{"GRASS":>8}{"grass>tree":>12}'
              f'{"BNS dens":>10}{"BNS cmass":>11}{"C3G cmass":>11}{"BNS anpp":>10}{"C3G anpp":>10}'
              f'{"cells BNS>.01":>14}{">.05":>7}')
        for y in sorted(F):
            f = F[y]; w = f['w']
            gt = np.mean(np.array(f['GRASSFPC#']) > np.array(f['TREEFPC#'])) * 100
            bns = np.array(f['BNS#'])
            g = lambda src, c: (src[y][c] / src[y]['w']) if (y in src and c in src[y]) else float('nan')
            print(f'  {y:6d}{f["BNS"]/w:9.4f}{f["C3G"]/w:9.4f}{f["TREEFPC"]/w:8.4f}{f["GRASSFPC"]/w:8.4f}'
                  f'{gt:11.1f}%{g(D,"BNS"):10.5f}{g(C,"BNS"):11.4f}{g(C,"C3G"):11.4f}'
                  f'{g(A,"BNS"):10.5f}{g(A,"C3G"):10.5f}{int((bns>0.01).sum()):14d}{int((bns>0.05).sum()):7d}')
        print()
