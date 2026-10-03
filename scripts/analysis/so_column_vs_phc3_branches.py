"""The Southern Ocean water column against PHC3: stratification and the deep heat reservoir.

Area-weighted mean profiles of potential temperature and salinity over open-ocean nodes
(no cavity) in 60-78S, for PHC3 (on the same 220509-node mesh), CONTROL over Y0-Y1 and
each branch, annual and JJA, whole band and three sectors (Weddell 300E-20E, Indian
20-150E, Pacific 150-300E).  Two stratification numbers per profile: the density step
between 0-50 m and 200-500 m (sigma-0 from a linear EOS, enough for a difference), split
into its temperature and salinity parts.  A model with too small a step convects; the
split says whether the halocline is too weak (salt) or the thermocline (heat).

Also the heat content anomaly against PHC3 in 100-500 m and 500-1500 m, in K, as the
reservoir the winter convection vents.

Usage:  CONTROL=PI200 RUNS=PI200_spp,PI200_h0 Y0=1565 Y1=1579 \
        python3 scripts/analysis/so_column_vs_phc3_branches.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
CONTROL = os.environ.get('CONTROL', 'PI200')
RUNS = os.environ.get('RUNS', 'PI200_spp,PI200_h0').split(',')
Y0, Y1 = int(os.environ.get('Y0', 1565)), int(os.environ.get('Y1', 1579))
MESH = '/work/ab0246/a270092/input/fesom2/core3/fesom.mesh.diag.nc'
PHC = '/work/ab0246/a270092/postprocessing/climatologies/CORE3_220509'
ALPHA, BETA = 0.6e-4, 0.78                                   # d(sigma)/dT [kg/m3/K] near 0 C, d(sigma)/dS [kg/m3/psu]

with xr.open_dataset(MESH) as m:
    lon = m['lon'].values % 360; lat = m['lat'].values
    area = m['nod_area'].values[0]; ul = m['ulevels_nod2D'].values; nlv = m['nlevels_nod2D'].values
    z = m['nz1'].values                                    # layer midpoints (47)
ok = (lat >= -78) & (lat <= -60) & (ul == 1)
SECT = {'all': ok, 'Weddell 300-20E': ok & ((lon >= 300) | (lon < 20)), 'Indian 20-150E': ok & (lon >= 20) & (lon < 150),
        'Pacific 150-300E': ok & (lon >= 150) & (lon < 300)}


def prof(a, k):
    """area-weighted mean over nodes k, per level, ignoring levels below the bottom (nan)"""
    w = area[k][:, None] * np.isfinite(a[k])
    return np.nansum(a[k] * w, 0) / np.where(w.sum(0) > 0, w.sum(0), np.nan)


def load_phc(v):
    with xr.open_dataset(f'{PHC}/{v}.fesom.1958.nc', decode_times=False) as d:
        a = np.squeeze(d[v].values)[:, :47].astype('f8')
    a[np.abs(a) > 1e30] = np.nan
    lev = np.arange(47)[None, :]
    return np.where(lev < (nlv - 1)[:, None], a, np.nan)


def load_run(arm, v, months):
    acc = 0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{arm}/outdata/fesom/{v}.fesom.{y}.nc', decode_times=False) as d:
            a = d[v].isel(time=months).mean('time').values.astype('f8')
        acc = acc + a / (Y1 - Y0 + 1)
    acc[np.abs(acc) > 1e30] = np.nan
    lev = np.arange(47)[None, :]
    return np.where(lev < (nlv - 1)[:, None], acc, np.nan)


P = {'T': load_phc('temp'), 'S': load_phc('salt')}
runs = {CONTROL: None}
for a in RUNS: runs[a] = None
DATA = {}
for arm in runs:
    DATA[arm] = {'ANN': {'T': load_run(arm, 'temp', list(range(12))), 'S': load_run(arm, 'salt', list(range(12)))},
                 'JJA': {'T': load_run(arm, 'temp', [5, 6, 7]), 'S': load_run(arm, 'salt', [5, 6, 7])}}
kz = lambda a, b: (z >= a) & (z < b)


def strat(T, S, k):
    t, s = prof(T, k), prof(S, k)
    up, lo = kz(0, 50), kz(200, 500)
    dT = np.nanmean(t[lo]) - np.nanmean(t[up]); dS = np.nanmean(s[lo]) - np.nanmean(s[up])
    return -ALPHA * dT * 1000, BETA * dS, dT, dS                  # sigma step from T, from S [kg/m3], and raw steps


print(f'60-78S open ocean, {Y0}-{Y1}.  Density step 200-500 m minus 0-50 m [kg/m3] split into T and S parts (positive = stable);')
print('layer temperatures [degC]; model minus PHC3 in brackets for the runs')
for sec, k in SECT.items():
    print(f'\n== {sec}  ({area[k].sum() / 1e12:.2f} M km2)')
    print(f'  {"":<14}{"dsig T":>8}{"dsig S":>8}{"total":>8}{"dS 0-50 vs 200-500":>20}{"T 0-50":>9}{"T100-500":>10}{"T500-1500":>11}{"S 0-50":>9}')
    for name, T, S in [('PHC3', P['T'], P['S'])] + [(f'{a} {s}', DATA[a][s]['T'], DATA[a][s]['S']) for a in runs for s in ('ANN', 'JJA')]:
        sT, sS, dT, dS = strat(T, S, k); t, s_ = prof(T, k), prof(S, k)
        lay = lambda f, a, b: np.nanmean(f[kz(a, b)])
        ref = (lambda f, a, b: '') if name == 'PHC3' else (lambda f, a, b: f' ({lay(f, a, b) - lay(prof(P["T"] if f is t else P["S"], k), a, b):+.2f})')
        print(f'  {name:<14}{sT:8.3f}{sS:8.3f}{sT + sS:8.3f}{dS:20.3f}{lay(t, 0, 50):9.2f}{ref(t, 0, 50):<8}{lay(t, 100, 500):7.2f}{ref(t, 100, 500):<8}'
              f'{lay(t, 500, 1500):7.2f}{ref(t, 500, 1500):<8}{lay(s_, 0, 50):7.3f}{ref(s_, 0, 50)}')

fig, axs = plt.subplots(2, 4, figsize=(15, 8), constrained_layout=True)
for j, (sec, k) in enumerate(SECT.items()):
    for i, v in enumerate(('T', 'S')):
        ax = axs[i, j]
        ax.plot(prof(P[v], k), z, 'k', lw=2, label='PHC3')
        for arm, ls in zip(runs, ('-', '--', ':', '-.')):
            ax.plot(prof(DATA[arm]['ANN'][v], k), z, ls, label=f'{arm} ANN')
            ax.plot(prof(DATA[arm]['JJA'][v], k), z, ls, lw=0.8, alpha=0.6, label=f'{arm} JJA')
        ax.set_ylim(1500, 0); ax.set_title(f'({chr(97 + 4 * i + j)}) {sec}: {"potential temperature [C]" if v == "T" else "salinity [psu]"}', fontsize=9)
        ax.grid(alpha=0.3)
        if v == 'S': ax.set_xlim(33.6, 34.9)
axs[0, 0].legend(fontsize=7)
fig.suptitle(f'60-78S water column, {Y0}-{Y1}, against PHC3')
out = os.path.join(REPO, 'plots', f'so_column_vs_phc3_{Y0}-{Y1}.png'); fig.savefig(out, dpi=110); print('saved', out)
