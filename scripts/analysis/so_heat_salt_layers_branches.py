"""Southern Ocean heat and salt by depth layer, branches against control, year by year.

Question: a lever that shuts off winter convection under the Antarctic pack (spp) should
leave the deep heat it used to vent trapped below a new halocline.  If so, the layers below
about 100 m warm year on year while the surface layer cools and freshens, and the stored
heat is a later risk (a convective release) rather than a fix.

Per run and year: annual-mean potential temperature and salinity, volume-weighted over
ocean cells in 60-78S, in the layers 0-100, 100-500, 500-1500 and 1500 m-bottom, from
FESOM's regular 0.5-degree output (temp.fesom.gr, salt.fesom.gr; nz is the layer mid-depth,
layer thickness from the midpoints).  Then the difference against CONTROL, per year and
for the second decade, with the heat change of each layer in 1e21 J.

Usage:  CONTROL=PI200 RUNS=PI200_spp,PI200_mle,PI200_h0,PI200_all3 Y0=1560 Y1=1579 \
        python3 scripts/analysis/so_heat_salt_layers_branches.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
CONTROL = os.environ.get('CONTROL', 'PI200')
RUNS = os.environ.get('RUNS', 'PI200_spp,PI200_mle,PI200_h0,PI200_all3').split(',')
Y0, Y1 = int(os.environ.get('Y0', 1560)), int(os.environ.get('Y1', 1579))
LAYERS = [('0-100', 0, 100), ('100-500', 100, 500), ('500-1500', 500, 1500), ('1500-bot', 1500, 1e5)]
RHO_CP = 1025.0 * 3990.0


def load(arm, v, y):
    with xr.open_dataset(f'{R}/{arm}/outdata/fesom/{v}.fesom.gr.{y}.nc', decode_times=False) as d:
        lat = d['lat'].values
        j = np.where((lat >= -78) & (lat <= -60))[0]
        a = d[v].isel(lat=slice(j[0], j[-1] + 1)).mean('time').values.astype('f8')   # (lat, lon, nz)
        z = d['nz'].values.astype('f8'); la = lat[j]
    a = np.where(np.abs(a) > 1e30, np.nan, a)
    return a, la, z


def layer_means(a, la, z):
    edges = np.concatenate([[0.0], 0.5 * (z[1:] + z[:-1]), [z[-1] + 0.5 * (z[-1] - z[-2])]])
    dz = np.diff(edges)
    area = np.cos(np.deg2rad(la))[:, None] * (0.5 * 111.195e3) ** 2 * np.ones(a.shape[1])[None, :]
    out = {}
    for n, top, bot in LAYERS:
        k = (z >= top) & (z < bot)
        vol = area[:, :, None] * dz[None, None, k] * np.isfinite(a[:, :, k])
        out[n] = (np.nansum(a[:, :, k] * vol) / vol.sum(), vol.sum())
    return out


D = {}
for arm in [CONTROL] + RUNS:
    D[arm] = {}
    for y in range(Y0, Y1 + 1):
        try:
            t, la, z = load(arm, 'temp', y); s, _, _ = load(arm, 'salt', y)
        except (FileNotFoundError, OSError):
            continue
        D[arm][y] = {'T': layer_means(t, la, z), 'S': layer_means(s, la, z)}
    print(f'{arm}: {len(D[arm])} years')

yrs = sorted(set.intersection(*[set(D[a]) for a in D]))
names = [n for n, _, _ in LAYERS]
print(f'\n60-78S ocean.  Control layer means, then branch minus control per year: T [K] / S [psu]')
for q, unit in (('T', 'K'), ('S', 'psu')):
    print(f'\n{q} [{unit}]   control ({CONTROL}) means ' + '  '.join(f'{n} {np.mean([D[CONTROL][y][q][n][0] for y in yrs]):.3f}' for n in names))
    print(f'  {"year":<6}' + ''.join(f'{a.replace(CONTROL + "_", ""):>34}' for a in RUNS))
    print(f'  {"":<6}' + ''.join(''.join(f'{n:>8}' for n in names) + '  ' for _ in RUNS))
    for y in yrs:
        row = f'  {y:<6}'
        for a in RUNS:
            row += ''.join(f'{D[a][y][q][n][0] - D[CONTROL][y][q][n][0]:+8.3f}' for n in names) + '  '
        print(row)

sec = [y for y in yrs if y >= yrs[0] + 10]
print(f'\nsecond decade {sec[0]}-{sec[-1]}: heat change against control [1e21 J] and T change [K], per layer')
for a in RUNS:
    cells = []
    for n in names:
        dT = np.mean([D[a][y]['T'][n][0] - D[CONTROL][y]['T'][n][0] for y in sec])
        vol = D[CONTROL][sec[0]]['T'][n][1]
        cells.append(f'{n} {dT * vol * RHO_CP / 1e21:+7.2f} ({dT:+.3f})')
    print(f'  {a:<12} ' + '   '.join(cells))
print('\ntrend of the branch-minus-control T difference over all common years [K per decade]')
for a in RUNS:
    tr = [np.polyfit(yrs, [D[a][y]['T'][n][0] - D[CONTROL][y]['T'][n][0] for y in yrs], 1)[0] * 10 for n in names]
    print(f'  {a:<12} ' + '   '.join(f'{n} {t:+.3f}' for n, t in zip(names, tr)))
