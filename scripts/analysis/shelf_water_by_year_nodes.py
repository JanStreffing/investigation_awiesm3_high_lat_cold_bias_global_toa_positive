"""Antarctic shelf water year by year: is the meltwater cut making the shelf saltier and denser?

Open-ocean nodes south of 65S with bottom shallower than DMAX (default 1000 m), no cavity,
area-weighted: annual 0-50 m and bottom-layer salinity and temperature, JJA mixed layer
(MLD2), and the fraction of shelf area where the annual-mean bottom water is denser than
sigma-0 27.8 (a dense-shelf-water proxy, linear EOS).  One line per year, so a run pair
can be compared while a leg is still running.

Usage:  ROOT=... ARM=PI200_view1600 Y0=1600 Y1=1607 MESH=<mesh diag> python3 scripts/analysis/shelf_water_by_year_nodes.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ['ROOT']; ARM = os.environ['ARM']; MESH = os.environ['MESH']; Y0, Y1 = int(os.environ['Y0']), int(os.environ['Y1']); DMAX = float(os.environ.get('DMAX', 1000))
with xr.open_dataset(MESH) as m:
    lat = m['lat'].values; nlv = m['nlevels_nod2D'].values; z = m['nz1'].values; zb = np.abs(m['zbar_n_bottom'].values)
    ul = m['ulevels_nod2D'].values if 'ulevels_nod2D' in m else np.ones(len(lat), int); area = m['nod_area'].values[0]
k = (lat < -65) & (ul == 1) & (zb < DMAX) & (zb > 50)
w = area[k]; ku = (z >= 0) & (z < 50)
print(f'{ARM}: shelf nodes south of 65S, bottom 50-{DMAX:.0f} m: {k.sum()} nodes, {w.sum() / 1e12:.2f} M km2')
print(f'  {"year":<6}{"S 0-50":>8}{"T 0-50":>8}{"S bottom":>10}{"T bottom":>10}{"sig0 bot":>10}{"dense frac":>11}{"MLD JJA":>9}')
for y in range(Y0, Y1 + 1):
    try:
        with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/temp.fesom.{y}.nc', decode_times=False) as d: T = d['temp'].mean('time').values.astype('f8')
        with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/salt.fesom.{y}.nc', decode_times=False) as d: S = d['salt'].mean('time').values.astype('f8')
        with xr.open_dataset(f'{R}/{ARM}/outdata/fesom/MLD2.fesom.{y}.nc', decode_times=False) as d: M = np.abs(d['MLD2'].isel(time=[5, 6, 7]).mean('time').values.astype('f8'))
    except (FileNotFoundError, OSError):
        continue
    if T.shape[0] != len(lat): T = T.T; S = S.T
    if T.shape[0] == 0 or not np.isfinite(T[k]).any(): continue
    T[np.abs(T) > 100] = np.nan; S[np.abs(S) > 100] = np.nan
    idx = np.arange(len(lat))[k]; ib = np.clip(nlv[k] - 2, 0, len(z) - 1)
    Tb = T[idx, ib]; Sb = S[idx, ib]
    Tu = np.nanmean(T[k][:, ku], 1); Su = np.nanmean(S[k][:, ku], 1)
    sig = 1000 * (0.6e-4 * (-Tb) * 0 + 0) + (-0.06 * Tb + 0.78 * Sb + 1000 * 0)     # linear sigma-0 [kg/m3] offsets cancel in the threshold below
    sig0 = 0.78 * (Sb - 34.5) - 0.06 * (Tb + 1.0) + 27.83                              # linear about (34.5 psu, -1 C): sigma-0 27.83
    av = lambda f: float(np.nansum(f * w * np.isfinite(f)) / np.sum(w * np.isfinite(f)))
    dense = float(np.sum(w[sig0 >= 27.8]) / w.sum())
    print(f'  {y:<6}{av(Su):8.3f}{av(Tu):8.2f}{av(Sb):10.3f}{av(Tb):10.2f}{av(sig0):10.3f}{dense:11.2f}{av(M[k]):9.0f}')
