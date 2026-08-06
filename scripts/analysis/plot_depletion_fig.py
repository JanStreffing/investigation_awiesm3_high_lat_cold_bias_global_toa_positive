"""Figure for the report: why the tanh broke the winter, and what replaces it.

Panel (a)  January snow-cover DISTRIBUTION across the Siberian box, N1 (as-released,
           clipping) against N2 (tanh).  The two have the SAME box-mean cover to within
           0.002 and the same total exposed area -- but N1 puts 93.5% of the box at
           cover exactly 1.0 while N2 dumps a quarter of it into 0.99-0.999.  That is
           the whole mechanism: a sliver of bare ground in every cell instead of a few
           genuinely bare cells, and each sliver couples that cell's soil to the air.

Panel (b)  The depletion curve itself, against 36492 RIHMI snow-course surveys, for a
           fresh-snow and an aged-snow density class.  tanh cannot reach 1 at any
           parameters; the as-released ramp reaches 1 but ignores density; the fitted
           form does both.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from runs import RT, LSMF

OBSD = '/work/ab0246/a270092/obs/RIHMI-WDC/data'
OUT = ('/work/ab0246/a270092/postprocessing/'
       'investigation_awiesm3_high_lat_cold_bias_global_toa_positive/report/plots/'
       'snow_depletion_mechanism.png')
YEARS = range(1876, 1896)
DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
DOY_MONTH = np.repeat(np.arange(12), DPM)
FIT = {'field': (0.014, 4.70, 1.46), 'forest': (0.026, 4.70, 0.40)}   # d_c(200), MD, b


def cov_tanh(d, rho, z0=0.016, rn=100.0, m=1.6):
    return np.tanh(d / np.maximum(2.5 * z0 * (rho / rn) ** m, 1e-6))


def cov_rel(d):
    return np.clip(10.0 * d, 0, 1)


def cov_fit(d, rho, p):
    dc0, md, b = p
    return np.clip((d / np.maximum(dc0 * (rho / 200.0) ** md, 1e-9)) ** b, 0, 1)


# ---------------------------------------------------------------- panel (a) ---
with xr.open_dataset(LSMF) as ds:
    lsm = ds['lsm'].isel(time_counter=0).values
    lat, lon = ds['lat'].values, ds['lon'].values
LA = np.broadcast_to(lat[:, None], lsm.shape)
LO = np.broadcast_to(lon[None, :], lsm.shape)
sel = np.flatnonzero(((LA >= 55) & (LA <= 75) & (LO >= 60) & (LO <= 180) & (lsm > 0.5)).ravel())
w = np.cos(np.deg2rad(LA)).ravel()[sel]; W = w / w.sum()


def load(run, var, y):
    f = f'{RT}/{run}/outdata/oifs/atm_remapped_1d_{var}_1d_{y}-{y}.nc'
    if not os.path.exists(f):
        return None
    with xr.open_dataset(f, decode_times=False) as d:
        a = d[var].values
    if a.shape[0] == 366:
        a = np.delete(a, 59, axis=0)
    return a.reshape(a.shape[0], -1)[:, sel] if a.shape[0] == 365 else None


EDGES = np.array([0.0, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0001])
HIST, MEANC = {}, {}
for lab, run, p in (('N1  as-released (clips)', 'amip_N1_snowdiag', None),
                    ('N2  tanh (asymptotic)', 'amip_N2_snowdiag_scf', 'tanh')):
    h, mc, n = np.zeros(len(EDGES) - 1), 0.0, 0
    for y in YEARS:
        sd, rsn = load(run, 'sd', y), load(run, 'rsn', y)
        if sd is None or rsn is None:
            continue
        jan = DOY_MONTH == 0
        rho = np.maximum(rsn[jan], 1e-6); d = sd[jan] * 1000.0 / rho
        c = cov_tanh(d, rho) if p == 'tanh' else cov_rel(d)
        for i in range(len(EDGES) - 1):
            h[i] += np.average(((c >= EDGES[i]) & (c < EDGES[i + 1])).mean(axis=0), weights=W)
        mc += np.average(c.mean(axis=0), weights=W); n += 1
    HIST[lab], MEANC[lab] = h / n, mc / n
    print(f'{lab}: mean cover {mc/n:.4f}  f_full {h[-1]/n:.4f}')

# ---------------------------------------------------------------- panel (b) ---
# Course type 0 == RIHMI type 1 == pole/field, which is what the FIELD parameters were
# fitted to.  The first version of this figure plotted field parameters against BOTH
# course types pooled, which made the curve look wrong when it was the comparison that
# was wrong.  Compare like with like.
with xr.open_dataset(f'{OBSD}/RIHMI-WDC_snmar.nc') as ds:
    la, lo = ds['lat'].values, ds['lon'].values
    k = (la >= 55) & (la <= 75) & (lo >= 60) & (lo <= 180)
    oc = ds['fraction_of_the_snow_course_covered_by_snow'].isel(station=k).values[:, :, 0] / 10.0
    od = ds['snow_depth_mean'].isel(station=k).values[:, :, 0] / 100.0
    orho = ds['snow_density'].isel(station=k).values[:, :, 0] * 1000.0

# ------------------------------------------------------------------- draw -----
fig, ax = plt.subplots(1, 2, figsize=(13.0, 4.6))
lbl = ['0–0.5', '0.5–0.8', '0.8–0.9', '0.9–0.95', '0.95–0.99', '0.99–0.999', '= 1.000']
x = np.arange(len(lbl)); wdt = 0.38
for i, (lab, col) in enumerate(zip(HIST, ('#2c7fb8', '#d95f0e'))):
    ax[0].bar(x + (i - 0.5) * wdt, HIST[lab], wdt, label=f'{lab}\nmean cover {MEANC[lab]:.3f}',
              color=col, edgecolor='k', linewidth=0.4)
ax[0].set_xticks(x); ax[0].set_xticklabels(lbl, rotation=30, ha='right', fontsize=8)
ax[0].set_ylabel('area fraction of the Siberian box')
ax[0].set_title('(a) January snow cover: identical mean, opposite state', fontsize=10)
ax[0].legend(fontsize=7.5, loc='upper left'); ax[0].grid(axis='y', alpha=0.3)
ax[0].annotate('every cell slightly bare\n$\\rightarrow$ no soil insulated',
               xy=(5, HIST['N2  tanh (asymptotic)'][5]), xytext=(2.4, 0.55), fontsize=8,
               arrowprops=dict(arrowstyle='->', lw=0.9))

dd = np.linspace(0.005, 0.6, 300)
for j, (rho0, ttl, ls) in enumerate(((130.0, 'fresh snow, $\\rho=130$', '-'),
                                     (285.0, 'aged spring snow, $\\rho=285$', '--'))):
    m = np.isfinite(oc) & np.isfinite(od) & np.isfinite(orho) & (np.abs(orho - rho0) < 30)
    be = np.array([.03, .06, .09, .13, .18, .25, .35, .5, .8])
    bx, by = [], []
    for i in range(len(be) - 1):
        s = m & (od >= be[i]) & (od < be[i + 1])
        if s.sum() >= 25:
            bx.append(od[s].mean()); by.append(oc[s].mean())
    ax[1].plot(bx, by, 'ko', ms=6, mfc='none' if j else 'k', mew=1.4,
               label=f'RIHMI field courses, {ttl}')
    ax[1].plot(dd, cov_fit(dd, rho0, FIT['field']), color='#238b45', ls=ls, lw=2,
               label='fitted  min$(1,(d/d_c)^b)$' if j == 0 else None)
    ax[1].plot(dd, cov_tanh(dd, rho0), color='#d95f0e', ls=ls, lw=1.6,
               label='tanh (cannot reach 1)' if j == 0 else None)
ax[1].plot(dd, cov_rel(dd), color='#2c7fb8', lw=1.6, label='as-released (no $\\rho$ term)')
ax[1].axhline(1.0, color='gray', lw=0.7, ls=':')
ax[1].set_xlim(0, 0.6); ax[1].set_ylim(0.3, 1.04)
ax[1].set_xlabel('snow depth $d$ [m]'); ax[1].set_ylabel('snow cover fraction')
ax[1].set_title('(b) depletion curve vs RIHMI field courses (solid $\\rho$=130, open $\\rho$=285)',
                fontsize=10)
ax[1].legend(fontsize=7.2, loc='lower right'); ax[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUT, dpi=155)
print('wrote', OUT)
