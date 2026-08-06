"""RIHMI-WDC station validation of the three quantities the snow scheme controls.

WHY THIS EXISTS.  The N/O series showed that Siberian winter soil temperature falls
17-25 K when ECE_SNOW_SCF is on, while DJF snow cover, SWE and depth are unchanged
to within 0.002 / 1 kg m-2 / 0.06 m.  Every box-mean explanation has failed.  The
surviving hypothesis is that the as-released ramp min(1,10d) CLIPS -- so most cells
sit at cover exactly 1.0 and their soil is perfectly insulated -- whereas tanh(d/L)
is asymptotic and can never report 1.0 anywhere, so the same box-mean cover is built
from every cell being slightly bare, and every cell's soil is then dragged toward air
temperature through the exposed tiles at lambda_sk = 10 W m-2 K-1.

That hypothesis makes an OBSERVATIONAL claim, not just a numerical one: real Siberian
midwinter snow cover reaches 100%, it does not asymptote to 96%.  RIHMI-WDC can test
it directly -- snow_cover_degree is the station analogue of ZCVS.

Three questions, three datasets:
  1. snow.nc   -- is observed DJF cover 10/10?  What fraction of station-days are at
                  complete cover?  This is the observational counterpart of f_full and
                  it decides whether a saturation cut is physics or a fudge.
  2. tpg.nc    -- observed soil temperature under NATURAL COVER, 1963-2024, 110
                  stations, 12 depths.  N1 gives 265.0 K and N2 gives 240.4 K in
                  January; only one of those can be real.  Shallow sensors (2-15 cm)
                  have no winter data, so 0.2 m and deeper are used, matched against
                  stl2 (7-28 cm, centre ~0.175 m).
  3. snmar.nc  -- snow-course density.  The tanh length scale is
                  L = 2.5*z0*(rho/rho_new)^m, so a density bias is a cover bias.  The
                  model carries 170-200 kg m-3 in midwinter; if observations say 220-260
                  the scheme is being driven off a biased density in the first place.

CAVEAT ON PERIOD.  The model runs are pre-industrial (1872-1915); the observations are
1963-2024.  That mismatch is worth a few K on air temperature and is stated with every
number below -- it is not worth the 25 K the winter soil discrepancy has to explain,
so the comparison is still decisive for THAT question, and only that one.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

OBSD = '/work/ab0246/a270092/obs/RIHMI-WDC/data'
RUNS = [('N1 ref off', 'amip_N1_snowdiag',     None),
        ('N2 cur',     'amip_N2_snowdiag_scf', (0.016, 100.0, 1.6)),
        ('O1 m4',      'amip_O1_scf_m4',       (0.018, 170.0, 4.0)),
        ('O2 m3',      'amip_O2_scf_m3',       (0.018, 170.0, 3.0))]
YEARS = range(1876, 1896)
DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
DOY_MONTH = np.repeat(np.arange(12), DPM)
BOX = dict(lat0=55, lat1=75, lon0=60, lon1=180)
RQSNCR_INV = 10.0


def cover(depth, rho, p):
    if p is None:
        return np.clip(RQSNCR_INV * depth, 0.0, 1.0)
    z0, rho_new, m = p
    scale = 2.5 * z0 * (np.maximum(rho, 50.0) / rho_new) ** m
    c = np.tanh(depth / np.maximum(scale, 1e-6))
    return np.where((depth > 1e-6) & (depth * rho > 1e-6), np.clip(c, 0, 1), 0.0)


def in_box(la, lo):
    return (la >= BOX['lat0']) & (la <= BOX['lat1']) & (lo >= BOX['lon0']) & (lo <= BOX['lon1'])


def tenths(a):
    """RIHMI reports cover in tenths, 0..10, with 99 as the MISSING sentinel.

    snow.nc carries 616663 values of 99 out of 5.4M.  Scaling those to 9.9 instead
    of masking them inflates the DJF mean to 1.97 on a quantity bounded by 1 and
    makes a quarter of Siberian midwinter look bare.  Anything above 10 is missing.
    """
    a = np.where(a > 10.0, np.nan, a)
    return a / 10.0


print(__doc__.split('CAVEAT ON PERIOD')[0])
print('=' * 94)

# ================================================================= 1. COVER ====
with xr.open_dataset(f'{OBSD}/RIHMI-WDC_snow.nc') as d:
    la, lo = d['lat'].values, d['lon'].values
    k = in_box(la, lo)
    cov = tenths(d['snow_cover_degree'].isel(station=k).values)   # (st, time)
    sdo = d['snow_depth'].isel(station=k).values                  # cm
    tm = d['time'].values
mo = tm.astype('datetime64[M]').astype(int) % 12
print(f'\n1. OBSERVED SNOW COVER DEGREE -- {k.sum()} stations in the Siberian box\n')
print(f'  {"month":8s}{"n obs":>10s}{"mean":>9s}{"median":>9s}{"f=10/10":>10s}{"f>=0.95":>10s}{"f<0.5":>9s}{"depth cm":>10s}')
obs_cov = {}
for m in (8, 9, 10, 11, 0, 1, 2, 3, 4):
    c = cov[:, mo == m].ravel(); c = c[np.isfinite(c)]
    s = sdo[:, mo == m].ravel(); s = s[np.isfinite(s)]
    if not c.size:
        continue
    obs_cov[m] = c.mean()
    print(f'  {MON[m]:8s}{c.size:10d}{c.mean():9.3f}{np.median(c):9.3f}'
          f'{(c >= 0.999).mean():10.3f}{(c >= 0.95).mean():10.3f}{(c < 0.5).mean():9.3f}'
          f'{(s.mean() if s.size else np.nan):10.1f}')
print('\n  f=10/10 is the observational counterpart of f_full: the fraction of')
print('  station-days at COMPLETE cover, which tanh(d/L) can never produce.')

# ============================================================ 2. SOIL TEMP ====
with xr.open_dataset(f'{OBSD}/RIHMI-WDC_tpg.nc') as d:
    la, lo = d['lat'].values, d['lon'].values
    k = in_box(la, lo)
    ts = d['tsoil'].isel(station=k).values            # (st, time, depth) degC
    dep = d['depth'].values
    tm2 = d['time'].values
    stla, stlo = la[k], lo[k]
mo2 = tm2.astype('datetime64[M]').astype(int) % 12
print(f'\n\n2. OBSERVED SOIL TEMPERATURE under natural cover -- {k.sum()} stations, 1963-2024\n')
print(f'  {"depth m":>9s}' + ''.join(f'{MON[m]:>9s}' for m in (9, 10, 11, 0, 1, 2, 5, 6)) + f'{"n DJF":>10s}')
for i, z in enumerate(dep):
    row, nd = [], 0
    for m in (9, 10, 11, 0, 1, 2, 5, 6):
        a = ts[:, mo2 == m, i].ravel(); a = a[np.isfinite(a)]
        row.append(a.mean() if a.size else np.nan)
        if m in (11, 0, 1):
            nd += a.size
    print(f'  {z:9.2f}' + ''.join(f'{v:9.1f}' for v in row) + f'{nd:10d}')

# ============================================================== 3. DENSITY ====
with xr.open_dataset(f'{OBSD}/RIHMI-WDC_snmar.nc') as d:
    la, lo = d['lat'].values, d['lon'].values
    k = in_box(la, lo)
    rho_o = d['snow_density'].isel(station=k).values * 1000.0     # g/cm3 -> kg/m3
    dep_o = d['snow_depth_mean'].isel(station=k).values           # cm
    frc_o = tenths(d['fraction_of_the_snow_course_covered_by_snow'].isel(station=k).values)
    tm3 = d['time'].values
mo3 = tm3.astype('datetime64[M]').astype(int) % 12
print(f'\n\n3. SNOW-COURSE DENSITY AND COVER -- {k.sum()} stations\n')
print(f'  {"month":8s}{"n":>8s}{"rho kg/m3":>12s}{"depth cm":>10s}{"course cover":>14s}')
obs_rho = {}
for m in (9, 10, 11, 0, 1, 2, 3, 4):
    r = rho_o[:, mo3 == m].ravel(); r = r[np.isfinite(r) & (rho_o[:, mo3 == m].ravel() > 0)]
    s = dep_o[:, mo3 == m].ravel(); s = s[np.isfinite(s)]
    f = frc_o[:, mo3 == m].ravel(); f = f[np.isfinite(f)]
    if not r.size:
        continue
    obs_rho[m] = r.mean()
    print(f'  {MON[m]:8s}{r.size:8d}{r.mean():12.1f}{(s.mean() if s.size else np.nan):10.1f}'
          f'{(f.mean() if f.size else np.nan):14.3f}')

# ================================================== 4. MODEL AT THE STATIONS ==
with xr.open_dataset(LSMF) as d:
    mlat, mlon = d['lat'].values, d['lon'].values
    lsm = d['lsm'].isel(time_counter=0).values
ny, nx = lsm.shape
iy = np.abs(mlat[None, :] - stla[:, None]).argmin(axis=1)
ix = np.abs(((mlon[None, :] - stlo[:, None] + 180) % 360) - 180).argmin(axis=1)
flat = np.unique(iy * nx + ix)


def load(run, var, y):
    for pat in (f'atm_remapped_1d_{var}_1d_{y}-{y}.nc', f'atm_remapped_1d_{var}_{y}-{y}.nc'):
        f = f'{RT}/{run}/outdata/oifs/{pat}'
        if os.path.exists(f):
            with xr.open_dataset(f, decode_times=False) as d:
                a = d[var].values
            if a.shape[0] == 366:
                a = np.delete(a, 59, axis=0)
            return a.reshape(a.shape[0], -1)[:, flat] if a.shape[0] == 365 else None
    return None


print(f'\n\n4. MODEL AT THE SAME {flat.size} GRID CELLS (nearest cell to each soil station)\n')
print(f'  {"run":12s}{"var":>8s}' + ''.join(f'{MON[m]:>9s}' for m in (9, 10, 11, 0, 1, 2)))
M = {}
for lab, run, p in RUNS:
    acc = {k: np.zeros((12,)) for k in ('cover', 'full', 'rho', 'stl1', 'stl2', 'depth')}
    n = 0
    for y in YEARS:
        sd, rsn = load(run, 'sd', y), load(run, 'rsn', y)
        s1, s2 = load(run, 'stl1', y), load(run, 'stl2', y)
        if sd is None or rsn is None or s1 is None or s2 is None:
            continue
        rho = np.maximum(rsn, 1e-6)
        dpt = sd * 1000.0 / rho
        c = cover(dpt, rho, p)
        for m in range(12):
            j = DOY_MONTH == m
            acc['cover'][m] += c[j].mean(); acc['full'][m] += (c[j] >= 0.999).mean()
            acc['rho'][m] += rho[j].mean(); acc['depth'][m] += dpt[j].mean()
            acc['stl1'][m] += s1[j].mean(); acc['stl2'][m] += s2[j].mean()
        n += 1
    if not n:
        continue
    M[lab] = {k: v / n for k, v in acc.items()}
    for v, fmt, off in (('cover', '9.3f', 0), ('full', '9.3f', 0), ('rho', '9.1f', 0),
                        ('depth', '9.3f', 0), ('stl1', '9.1f', -273.15), ('stl2', '9.1f', -273.15)):
        print(f'  {lab if v == "cover" else "":12s}{v:>8s}'
              + ''.join(format(M[lab][v][m] + off, fmt) for m in (9, 10, 11, 0, 1, 2)))

# ------------------------------------------------------------------ verdict ---
print('\n\n5. THE THREE COMPARISONS\n')
d20 = int(np.abs(dep - 0.2).argmin())
o_s = {m: np.nanmean(ts[:, mo2 == m, d20]) for m in (11, 0, 1)}
print(f'  soil at {dep[d20]:.2f} m, DJF observed : {np.mean(list(o_s.values())):+7.1f} degC')
for lab in M:
    print(f'    {lab:12s} stl2 DJF         : '
          f'{np.mean([M[lab]["stl2"][m] for m in (11, 0, 1)]) - 273.15:+7.1f} degC'
          f'   bias {np.mean([M[lab]["stl2"][m] for m in (11,0,1)]) - 273.15 - np.mean(list(o_s.values())):+7.1f}')
print(f'\n  DJF cover observed              : {np.mean([obs_cov[m] for m in (11,0,1) if m in obs_cov]):7.3f}')
for lab in M:
    print(f'    {lab:12s} cover / f_full   : {np.mean([M[lab]["cover"][m] for m in (11,0,1)]):7.3f}'
          f' / {np.mean([M[lab]["full"][m] for m in (11,0,1)]):7.3f}')
print(f'\n  DJF density observed            : {np.mean([obs_rho[m] for m in (11,0,1) if m in obs_rho]):7.1f} kg/m3')
for lab in M:
    print(f'    {lab:12s} rsn DJF          : {np.mean([M[lab]["rho"][m] for m in (11,0,1)]):7.1f}'
          f'   bias {np.mean([M[lab]["rho"][m] for m in (11,0,1)]) - np.mean([obs_rho[m] for m in (11,0,1) if m in obs_rho]):+7.1f}')
print("""
  READING IT.  If observed DJF cover is ~1.0 with most station-days at 10/10 while the
  scheme runs show f_full ~ 0, the scheme is wrong about a quantity that is directly
  measured, and the saturation cut is required by observation rather than convenience.
  If observed soil at 0.2 m is near N1 and nowhere near N2, the scheme's winter soil is
  a defect, not a correction of a compensating error.  And if observed density is well
  above the model's, the tanh is being driven off a biased L on top of everything else.""")
