"""Monthly area-weighted sea-ice fraction (a_ice) in the Labrador Sea along the tuning lineage.

Regions (node indices and FESOM dual-cell areas from mesh.nc, in lab_regions.npz):
  labsea    52-66N, 65-45W  (basin incl. Labrador shelf and West Greenland coast)
  interior  56-62N, 60-50W  (the convection box used by labrador_2x2.py)
Usage: python labrador_ice_lineage.py [var=MLD2] <regions.npz> <out.csv> <exp> <outdata_fesom_dir> <y0> <y1> [...]
var defaults to a_ice. For any other variable (e.g. MLD2) masked values are excluded,
not zeroed, and the absolute value is averaged.
Daily or monthly a_ice files both work; daily are averaged by calendar month.
Missing (masked) values count as zero ice.
"""
import sys, calendar, numpy as np
from netCDF4 import Dataset
argv = sys.argv[1:]
VAR = 'a_ice'
if argv[0].startswith('var='):
    VAR = argv.pop(0)[4:]
reg = np.load(argv[0]); out = open(argv[1], 'a')
args = argv[2:]
for i in range(0, len(args), 4):
    exp, d, y0, y1 = args[i], args[i+1], int(args[i+2]), int(args[i+3])
    for y in range(y0, y1 + 1):
        try:
            with Dataset(f'{d}/{VAR}.fesom.{y}.nc') as nc:
                v = nc.variables[VAR][:]; a = v.filled(np.nan) if hasattr(v, 'filled') else v
        except OSError:
            print(f'{exp} {y} missing', flush=True); continue
        nt = a.shape[0]
        if nt in (365, 366):
            ml = [calendar.monthrange(y, m)[1] for m in range(1, 13)] if nt == 366 else [31,28,31,30,31,30,31,31,30,31,30,31]
            edges = np.cumsum([0] + ml)
        elif nt == 12:
            edges = np.arange(13)
        else:
            print(f'{exp} {y} unexpected nt={nt}', flush=True); continue
        for r in ('labsea', 'interior'):
            idx, w = reg[r + '_idx'], reg[r + '_w']
            # a_ice output masks ice-free nodes as missing (value==0 masking in the
            # XIOS output of these runs); every mesh node in both boxes is ocean, so
            # missing means no ice.
            x = np.asarray(a[:, idx], dtype='f8')
            if VAR == 'a_ice':
                x = np.nan_to_num(x, nan=0.0); x[np.abs(x) > 1e10] = 0.0
                s = (x * w).sum(axis=1) / w.sum()
            else:
                x[np.abs(x) > 1e10] = np.nan; x = np.abs(x); ok = np.isfinite(x)
                s = np.where(ok, x, 0).dot(w) / (ok * w).sum(axis=1)
            for m in range(12):
                out.write(f'{exp},{y},{m+1},{r},{np.mean(s[edges[m]:edges[m+1]]):.5f}\n')
        out.flush(); print(f'{exp} {y} ok', flush=True)
