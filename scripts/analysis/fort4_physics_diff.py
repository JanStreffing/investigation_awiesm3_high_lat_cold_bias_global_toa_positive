"""Compare two OpenIFS fort.4 files group by group, as the model actually ran them.

Parses every &GROUP ... / block into key = value pairs (keys upper-cased, whitespace
normalised) and prints, per group, keys set in only one file and keys set to different
values.  Groups that are resolution- or run-control-driven by construction (grid,
timestep, dynamics, I/O, restart, date) are listed separately so the physics differences
stand out.

Usage:  python3 scripts/analysis/fort4_physics_diff.py <fort.4 A> <fort.4 B> [labelA labelB]
"""
import re, sys
RUNCTL = {'NAMRES', 'NAMCT0', 'NAMCT1', 'NAMDIM', 'NAMGEM', 'NAMDYN', 'NAMDYNA', 'NAMRIP', 'NAMPAR0', 'NAMPAR1',
          'NAMIO_SERV', 'NAMFPC', 'NAMFPD', 'NAMFPG', 'NAMFPIOS', 'NAMOPH', 'NAMARG', 'NAMCVA', 'NAMGRIB', 'NAMDDH',
          'NAMSTA', 'NAMVV1', 'NAMVV2', 'NAMFA', 'NAMMCC', 'NAMIOMI', 'NAMCHK', 'NAMSCC', 'NAMTRAJ', 'NAMCOM',
          'NAMPPC', 'NAMFPPHY', 'NAMFPDYH', 'NAMFPDYP', 'NAMFPDYV', 'NAMFPDYT', 'NAMFPDYS', 'NAMFPDYF', 'NAMFPEZO'}


def parse(path):
    groups = {}; cur = None
    for raw in open(path, errors='replace'):
        line = raw.split('!')[0].rstrip()
        m = re.match(r'\s*&(\w+)', line)
        if m:
            cur = m.group(1).upper(); groups.setdefault(cur, {}); line = line[m.end():]
        if cur is None:
            continue
        for kv in re.findall(r"([A-Za-z_][\w%]*(?:\([^)]*\))?)\s*=\s*('[^']*'|[^,/]+)", line):
            groups[cur][kv[0].upper().replace(' ', '')] = kv[1].strip().rstrip(',').strip()
        if re.search(r'^\s*/\s*$', line) or line.strip().endswith('/'):
            cur = None
    return groups


def num(v):
    try:
        return float(v.replace('d', 'e').replace('D', 'E'))
    except ValueError:
        return v.strip().upper()


a, b = parse(sys.argv[1]), parse(sys.argv[2])
la, lb = (sys.argv[3], sys.argv[4]) if len(sys.argv) > 4 else ('A', 'B')
phys, ctl = [], []
for g in sorted(set(a) | set(b)):
    A, B = a.get(g, {}), b.get(g, {})
    lines = []
    for k in sorted(set(A) | set(B)):
        va, vb = A.get(k), B.get(k)
        if va is None or vb is None or num(va) != num(vb):
            lines.append(f'    {k:34s} {la}: {va if va is not None else "(unset)":>18s}   {lb}: {vb if vb is not None else "(unset)":>18s}')
    if lines:
        (ctl if g in RUNCTL else phys).append((g, lines))
print(f'==== PHYSICS / TUNING groups that differ ({la} vs {lb})')
for g, lines in phys:
    print(f'&{g}'); print('\n'.join(lines))
print(f'\n==== run-control / resolution groups that differ (by construction): {", ".join(g for g, _ in ctl)}')
