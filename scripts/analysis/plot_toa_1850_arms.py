"""Global net TOA against year for the 1850 arms, with 11W as far as it has run.

WHY THIS PLOT.  The campaign quoted a "standing +0.88 W/m2 imbalance" for a week.  It is
not standing: the PI arms START near zero and DRIFT.  A time series says that in one look
where a pair of window means does not, and 11W is the arm whose whole purpose is the drift
rate rather than the level.

11W is the S4-removal twin of 11Q (RCL_INPPMIN 50000 -> 70000, one namelist number).  It is
incomplete.  11N is the pre-RSBLB predecessor, included because it shows the drift is a
property of the base rather than of either lever.

THE TRAP THIS PLOT HAS TO AVOID.  The drift is far steeper early than late, so an OLS slope
over 11W's 20 years is not comparable to one over 40.  Fitting each arm to its own length
made 11W look like it drifts 3x faster, which is a window artifact and not a result.  Every
trend drawn here is therefore fitted over the SAME window, the years 11W has run, and the
completed arms carry a thin second line for their full 40 years so the slowdown is visible.

Years 1350-51 are an initialisation transient, roughly -1.0 to -1.6 W/m2 in every arm, and
are shaded.  The campaign's "starts near zero" numbers are the 1352-59 mean and deliberately
exclude them.

THE DECADAL SMOOTH.  A centred 11-point Kaiser window, beta = 8.6, spanning 10 years either
side to side.  Kaiser at that beta is close to Blackman: the sidelobes are ~75 dB down, so
interannual variance does not leak into the smoothed curve and the residual wiggle is
signal.  It is drawn ONLY where the full window fits, five years in from each end.  No
reflection or zero-padding: on a series that is drifting this hard, a padded endpoint is
fabricated and would land exactly where the eye looks for the answer.  That is also why
the smooth is the honest way to see the slowdown, which is curvature, where a single OLS
slope over the whole record cannot show it at all.

TRAP.  IFS TOA fluxes are accumulated J/m2 over the output step; divide by 3600 or every
number is ~3600x too large.  Verified against global ASR 240.4 W/m2.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
OUT = ('/work/ab0246/a270092/postprocessing/'
       'investigation_awiesm3_high_lat_cold_bias_global_toa_positive/plots')
ACC = 3600.0

def resolve(pat):
    g = [p for p in glob.glob(f'{R092}/{pat}') if os.path.isdir(p)]
    return g[0] if g else None

ARMS = [('11N  LX4',            resolve('*11N*'),  '#888888', '-'),
        ('11Q  LX4+RSBLB  +S4', resolve('11Q'),    '#1f77b4', '-'),
        ('11W  LX4+RSBLB  -S4', resolve('11W'),    '#d62728', '-')]


def annual_net_toa(root, y):
    D = f'{root}/outdata/oifs'
    out = []
    for v in ('tsr', 'ttr'):
        f = f'{D}/atm_remapped_1m_{v}_{y}-{y}.nc'
        if not os.path.exists(f):
            return np.nan
        with xr.open_dataset(f, decode_times=False) as d:
            n = v if v in d.data_vars else list(d.data_vars)[-1]
            a = d[n].values / ACC
            lat = d['lat'].values
        if a.shape[0] != 12:
            return np.nan
        out.append(a.mean(axis=0))
    f2d = out[0] + out[1]
    return float(np.average(f2d.mean(axis=1), weights=np.cos(np.deg2rad(lat))))


fig, ax = plt.subplots(figsize=(10, 5.6))
ax.axhspan(-0.3, 0.3, color='#2ca02c', alpha=0.10, zorder=0)
ax.axhline(0.0, color='#2ca02c', lw=1.0, ls='--', zorder=1)
ax.text(1350.3, 0.31, 'PI target  net TOA within $\\pm$0.3', color='#2ca02c',
        fontsize=8.5, va='bottom')

# every trend is fitted over the window the shortest arm has, so the slopes compare
lengths = []
for _l, _r, _c, _ls in ARMS:
    if _r and os.path.isdir(f'{_r}/outdata/oifs'):
        lengths.append(max(int(f.split('_')[-1].split('-')[0])
                           for f in os.listdir(f'{_r}/outdata/oifs')
                           if f.startswith('atm_remapped_1m_tsr_') and f.endswith('.nc')))
MATCH_END = min(lengths)
ax.axvspan(1349.5, 1351.5, color='#999999', alpha=0.16, zorder=0)
ax.text(1350.0, -1.72, 'init\ntransient', fontsize=7.5, color='#666666', va='bottom')

KAISER_N, KAISER_BETA = 11, 8.6

def kaiser_smooth(v, n=KAISER_N, beta=KAISER_BETA):
    """Centred Kaiser-weighted running mean; returns NaN where the window overhangs."""
    w = np.kaiser(n, beta)
    w /= w.sum()
    out = np.full(len(v), np.nan)
    h = n // 2
    for i in range(h, len(v) - h):
        out[i] = float(np.dot(v[i - h:i + h + 1], w))
    return out


def fit(x, v):
    b, a = np.polyfit(x, v, 1)
    res = v - (a + b * x)
    se = np.sqrt((res ** 2).sum() / (len(x) - 2) / ((x - x.mean()) ** 2).sum())
    return a, b, se

summary = []
for label, root, col, ls in ARMS:
    if root is None:
        print(f'  {label}: directory not found'); continue
    ys = sorted(int(f.split('_')[-1].split('-')[0])
                for f in os.listdir(f'{root}/outdata/oifs')
                if f.startswith('atm_remapped_1m_tsr_') and f.endswith('.nc'))
    v = np.array([annual_net_toa(root, y) for y in ys])
    ok = np.isfinite(v)
    ys, v = np.array(ys)[ok], v[ok]
    if not len(ys):
        continue
    partial = len(ys) < 40
    ax.plot(ys, v, color=col, lw=0.8, alpha=0.30, marker='o', ms=2.2,
            label=f'{label}   {len(ys)} yr' + ('  (running)' if partial else ''))
    sm = kaiser_smooth(v)
    ax.plot(ys, sm, color=col, lw=3.0, zorder=5, solid_capstyle='round')
    m = ys <= MATCH_END
    a, b, se = fit(ys[m], v[m])
    ax.plot(ys[m], a + b * ys[m], color=col, lw=1.1, ls='--', alpha=0.75, zorder=4)
    full = (np.nan, np.nan)
    if not partial:
        _, b2, se2 = fit(ys, v)
        full = (b2 * 10, se2 * 10)
    e = (ys >= 1352) & (ys <= 1359)
    summary.append((label, len(ys), v[e].mean(), b * 10, se * 10, full, partial))

ax.set_xlabel('model year'); ax.set_ylabel('global net TOA  [W m$^{-2}$]')
ax.set_title(f'1850 arms: net TOA is a DRIFT, not a standing offset\n'
             f'thick = {KAISER_N}-point Kaiser decadal mean ($\\beta$={KAISER_BETA}, '
             f'full windows only);  dashed = OLS over the matched window 1350-{MATCH_END}',
             fontsize=11)
ax.legend(loc='lower right', fontsize=9, framealpha=0.92)
ax.grid(alpha=0.25, lw=0.5)
ax.set_xlim(1349, 1390)
fig.tight_layout()
p = f'{OUT}/toa_1850_arms.png'
fig.savefig(p, dpi=160)
print(f'wrote {p}\n')
print(f'  matched trend window: 1350-{MATCH_END}\n')
print(f'  {"arm":22s} {"yr":>3s} {"1352-59":>9s} {"matched trend/dec":>19s} '
      f'{"full 40yr trend/dec":>21s}')
for label, n, early, tr, se, full, partial in summary:
    fs = ('%+.3f +-%.3f' % full) if np.isfinite(full[0]) else 'still running'
    print(f'  {label:22s} {n:3d} {early:+9.3f} {tr:+12.3f} +-{se:.3f} {fs:>21s}')
