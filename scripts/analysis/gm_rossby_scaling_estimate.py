"""What GM coefficient would FESOM's scaling options give, by latitude band, on this mesh?

Reproduces oce_fer_gm.F90's node-wise scaling from the run's own annual stratification:
  c1   = max(K_GM_cmin, (1/pi) * integral N dz)        first baroclinic gravity wave speed
  rosb = min(c1 / max(|f|, 1e-6), 200 km)              Rossby radius
  Fermi cutoff (scaling_Rossby): 1 / (1 + exp(-(min(reso/rosb, 5) - 1.5) / 0.15))
  resolution scaling (order 2):  sqrt(2 * area) / 100 km
with N^2 from a linear equation of state on FESOM's regular output (temp, salt) and the
node spacing from the mesh diag, binned to the same 0.5-degree grid.  Prints the
area-weighted mean factor and the resulting K_GM per band for the current setting, for
Rossby-only, and for both, at K_GM_max.  An estimate: FESOM computes N from its full EOS
and per node, and the Fermi function is sharp, so band means can differ from this by tens
of percent; the year-1 fer_K output of a test run is the check.

Usage:  ARM=PI200 YEAR=1575 KGM=2500 python3 scripts/analysis/gm_rossby_scaling_estimate.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4'); ARM = os.environ.get('ARM', 'PI200')
Y = int(os.environ.get('YEAR', 1575)); KGM = float(os.environ.get('KGM', 2500))
MESH = '/work/ab0246/a270092/input/fesom2/core3/fesom.mesh.diag.nc'
G, RHO0, ALPHA, BETA = 9.81, 1025.0, 1.7e-4, 7.6e-4
with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/temp.fesom.gr.{Y}.nc', decode_times=False) as d:
    T = d['temp'].mean('time').values.astype('f8'); lat = d['lat'].values; lon = d['lon'].values; z = d['nz'].values
with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/salt.fesom.gr.{Y}.nc', decode_times=False) as d:
    S = d['salt'].mean('time').values.astype('f8')
T[np.abs(T) > 100] = np.nan; S[np.abs(S) > 100] = np.nan
sig = RHO0 * (-ALPHA * T + BETA * S)
n2 = np.full(sig.shape, np.nan); n2[..., :-1] = (G / RHO0) * (sig[..., 1:] - sig[..., :-1]) / (z[1:] - z[:-1])[None, None, :]   # z is depth, positive down
dz = np.diff(np.concatenate([[0.0], 0.5 * (z[1:] + z[:-1]), [z[-1] + 0.5 * (z[-1] - z[-2])]]))
c1 = np.nansum(np.sqrt(np.maximum(np.nan_to_num(n2), 0)) * dz[None, None, :], -1) / np.pi
ocean = np.isfinite(T[..., 0]); c1 = np.where(ocean, np.maximum(c1, 0.1), np.nan)      # K_GM_cmin = 0.1
f = 2 * 7.2921e-5 * np.sin(np.deg2rad(lat))[:, None] * np.ones(len(lon))[None, :]
rosb = np.minimum(c1 / np.maximum(np.abs(f), 1e-6), 2e5)
with xr.open_dataset(MESH) as m:
    nlat = m['lat'].values; nlon = m['lon'].values % 360; area = m['nod_area'].values[0]
iy = np.clip(((nlat + 90) / 0.5).astype(int), 0, len(lat) - 1); ix = np.clip((nlon / 0.5).astype(int), 0, len(lon) - 1)
num = np.zeros(T.shape[:2]); den = np.zeros_like(num)
np.add.at(num, (iy, ix), area); np.add.at(den, (iy, ix), 1.0)
A = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)               # mean node area per cell
# fill cells without a node from the zonal mean of the band (coarse cells at 0.5 deg are common)
for j in range(len(lat)):
    if np.isfinite(A[j]).any(): A[j] = np.where(np.isfinite(A[j]), A[j], np.nanmean(A[j]))
reso = np.sqrt(2 * A)                                                        # mesh_resolution proxy
fermi = 1.0 / (1.0 + np.exp(-(np.minimum(reso / rosb, 5.0) - 1.5) / 0.15))
resol = np.sqrt(2 * A) / 1e5
W = np.cos(np.deg2rad(lat))[:, None] * np.ones(len(lon))[None, :]
print(f'{ARM} {Y}: GM scaling by band (ocean cells, area-weighted), K_GM_max {KGM:.0f}')
print(f'  {"band":<10}{"c1 m/s":>8}{"Rossby km":>10}{"reso km":>8}{"resol fac":>10}{"Fermi fac":>10}{"K now":>7}{"K Rossby only":>14}{"K both":>8}')
for lo, hi in ((-78, -60), (-60, -45), (-45, -30), (-30, 30), (30, 45), (45, 60), (60, 80)):
    k = ocean & (lat[:, None] >= lo) & (lat[:, None] < hi) & np.isfinite(fermi) & np.isfinite(resol)
    av = lambda x: float(np.average(x[k], weights=W[k]))
    print(f'  {lo:>4}..{hi:<4}{av(c1):8.2f}{av(rosb) / 1e3:10.1f}{av(reso) / 1e3:8.1f}{av(resol):10.2f}{av(fermi):10.2f}'
          f'{KGM * av(resol):7.0f}{KGM * av(fermi):14.0f}{KGM * av(fermi * resol):8.0f}')
