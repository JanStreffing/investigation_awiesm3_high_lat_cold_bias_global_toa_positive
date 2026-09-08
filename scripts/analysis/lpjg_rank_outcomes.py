#!/usr/bin/env python3
"""Per-rank outcome table for an LPJ-GUESS -islpjgspinup run.

In spin-up mode every rank works through its cells sequentially and a fail() on one
cell kills the rank, so a state directory with N files says nothing about how many
cells were actually simulated.  This parses every run<r+1>/guess<r>.log and reports,
per cell: simulated to the end, died (with which message), or never reached.

Usage: lpjg_rank_outcomes.py <work_dir> [nranks] [--peat <peat_frac.txt>] [--csv out.csv]
"""
import os
import re
import sys
from collections import Counter
import numpy as np

work = sys.argv[1]
nranks = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 1280
peat_file = sys.argv[sys.argv.index('--peat') + 1] if '--peat' in sys.argv else None
csv_out = sys.argv[sys.argv.index('--csv') + 1] if '--csv' in sys.argv else None

grid = np.loadtxt(os.path.join(work, 'ece_gridlist_TCO95.txt'))
ncell = len(grid)
ifs_index = grid[:, 3].astype(int)          # 1-based index into the 40320-cell grid
by_ifs = {ifs_index[i]: i for i in range(ncell)}

re_read = re.compile(r'About to read T2M/tas for year=\d+, ifs_index=(\d+)')
re_saved = re.compile(r'has saved the state after spinup for cell with lon=([-\d.]+), lat=([-\d.]+)')
FAILS = [('LUH3: all 12 required state fields are missing', 'luh3_missing'),
         ('Cannot initialize soil for newly created', 'soil_init'),
         ('Fixed peat invariant', 'peat_invariant'),
         ('LPJG: configured restart could not be loaded', 'no_state')]

status = np.full(ncell, 'unreached', dtype=object)
rank_summary = []
for r in range(nranks):
    f = os.path.join(work, f'run{r + 1}', f'guess{r}.log')
    if not os.path.exists(f):
        rank_summary.append((r, 'nolog'))
        continue
    cur, died = None, None
    for line in open(f, errors='replace'):
        m = re_read.search(line)
        if m:
            cur = by_ifs.get(int(m.group(1)))
            continue
        if re_saved.search(line):
            if cur is not None:
                status[cur] = 'ok'
            continue
        for key, tag in FAILS:
            if line.startswith(key):
                died = tag
        if died:
            if cur is not None and status[cur] != 'ok':
                status[cur] = died
            break
    rank_summary.append((r, died or 'ok'))

print('cells:', ncell)
print('cell outcomes:', dict(Counter(status)))
print('rank outcomes:', dict(Counter(s[1] for s in rank_summary)))

if peat_file:
    peat = {}
    with open(peat_file) as fh:
        fh.readline()
        for line in fh:
            p = line.split()
            if len(p) >= 3:
                peat[(round(float(p[0]), 3), round(float(p[1]), 3))] = float(p[2])
    pf = np.array([peat.get((round(grid[i, 0], 3), round(grid[i, 1], 3)), np.nan) for i in range(ncell)])
    for s in ['ok', 'soil_init', 'luh3_missing', 'unreached']:
        sel = status == s
        if sel.any():
            v = pf[sel]
            print(f'{s:13s} n={sel.sum():5d}  peat>0: {np.nanmean(v > 0):.3f}  peat mean {np.nanmean(v):.4f}  no-map {np.isnan(v).sum()}')

if csv_out:
    with open(csv_out, 'w') as fh:
        fh.write('lon,lat,ifs_index,status\n')
        for i in range(ncell):
            fh.write(f'{grid[i,0]:.6f},{grid[i,1]:.6f},{ifs_index[i]},{status[i]}\n')
    print('wrote', csv_out)
