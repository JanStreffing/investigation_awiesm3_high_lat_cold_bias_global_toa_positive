"""Surface heat entry (-fh) into the Southern Ocean by latitude band, per year, for one run.
Splits the SO (< 34S) row of heat_pathways.py so the GM 1500 - 2500 difference can be placed
in the mode-water belt (40-55S) or further south. Area from fesom.mesh.diag.nc nod_area(1,:).
Usage: python so_sfc_entry_bands.py <mesh.diag.nc> <outdata_fesom> <exp> <y0> <y1> <out.csv>
Rows: exp,year,band,W
"""
import sys, numpy as np
from netCDF4 import Dataset
diag, d, exp, y0, y1, outf = sys.argv[1:7]; y0, y1 = int(y0), int(y1)
with Dataset(diag) as nc:
    A = np.asarray(nc.variables['nod_area'][0], 'f8'); lat = np.asarray(nc.variables['lat'][:], 'f8')
if np.abs(lat).max() < 1.6: lat = np.rad2deg(lat)
B = [(-40, -34), (-45, -40), (-50, -45), (-55, -50), (-60, -55), (-65, -60), (-90, -65)]
out = open(outf, 'a')
for y in range(y0, y1 + 1):
    with Dataset(f'{d}/fh.fesom.{y}.nc') as nc:
        v = nc.variables['fh']; acc = np.zeros(v.shape[1])
        for m in range(v.shape[0]):
            a = np.asarray(v[m], 'f8'); a[~np.isfinite(a) | (np.abs(a) > 1e5)] = 0.0; acc += a
        fh = acc / v.shape[0]
    for lo, hi in B:
        m = (lat >= lo) & (lat < hi)
        out.write(f'{exp},{y},{lo}..{hi},{-(fh[m] * A[m]).sum():.6e}\n')
    out.flush(); print(exp, y, 'ok', flush=True)
