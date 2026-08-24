"""Gregory regression on the 1850 arms: does PI reach net TOA zero at a plausible PI temperature?

THE QUESTION.  Not the drift rate, and not the plateau level on its own.  Net TOA and global
mean temperature are both still moving, so the thing that decides whether this configuration
is usable for a PI spin-up is where the two cross zero together: regress net TOA N on global
mean T2m, and the x-intercept is the temperature the model equilibrates at.  If that lands
near the PI value the configuration is sound; if it lands well above, the model will spin up
too warm no matter how long it runs.

WHAT A GREGORY PLOT NORMALLY IS, AND WHY THIS ONE IS DIFFERENT.  The classic construction
regresses N against dT after an ABRUPT forcing step (4xCO2), where the slope is the climate
feedback parameter and the intercept on the N axis is the forcing.  There is no step here.
These arms are relaxing from an imposed initial state at fixed 1850 forcing, so the slope is
an effective relaxation rate that mixes radiative feedback with ocean and sea-ice adjustment,
NOT a clean feedback parameter.  The x-intercept is still the equilibrium the run is heading
for, which is the quantity asked for.  The slope is reported but should not be quoted as
lambda.

A SIGN TO WATCH.  For a damped system N falls as T rises, giving a NEGATIVE slope and a
finite intercept.  These arms have been drifting the other way, N rising while T rises, which
is the ice-albedo feedback outrunning the damping.  A positive slope means no equilibrium in
the fitted window and the intercept is meaningless; the code refuses to print one in that
case rather than reporting a spurious crossing.

WINDOWS.  Years 1350-51 are an initialisation transient and are always excluded.  The fit is
reported over the full post-transient record and over a LATE window separately, because the
decadal Kaiser mean shows N flattening near +0.9 by about 1370: early and late are different
regimes and one line through both is not a relaxation.  The late window is the SAME calendar
years for every arm, not each arm's own second half.  Using each arm's own half compared
11W's approach phase against the completed arms' plateau and made 11W look like the only one
relaxing, which is a window artifact of exactly the kind that already had to be corrected out
of the TOA time-series plot.

WHEN AN INTERCEPT IS REFUSED.  A crossing is printed only if the slope is significantly
negative, |b| > 1.96*se.  A slope that is merely negative but consistent with zero puts the
crossing anywhere from just off the data to infinity, and an early version of this script
duly reported 18.52 C for 11N off a slope of -0.168 +- 0.422.  That number was meaningless
and is precisely the kind of extrapolation the reader would have believed.

TRAP.  IFS TOA fluxes are accumulated J/m2 over the output step; divide by 3600.
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
TRANSIENT_END = 1351

def resolve(p):
    g = [x for x in glob.glob(f'{R092}/{p}') if os.path.isdir(x)]
    return g[0] if g else None

ARMS = [('11N  LX4',            resolve('*11N*'), '#888888'),
        ('11Q  LX4+RSBLB  +S4', resolve('11Q'),   '#1f77b4'),
        ('11W  LX4+RSBLB  -S4', resolve('11W'),   '#d62728')]


def gm(root, var, y):
    f = f'{root}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(f):
        return None, None
    with xr.open_dataset(f, decode_times=False) as d:
        n = var if var in d.data_vars else list(d.data_vars)[-1]
        a = d[n].values
        lat = d['lat'].values
    if a.shape[0] != 12:
        return None, None
    return a, lat


def series(root):
    ys = sorted(int(f.split('_')[-1].split('-')[0])
                for f in os.listdir(f'{root}/outdata/oifs')
                if f.startswith('atm_remapped_1m_2t_') and f.endswith('.nc'))
    Y, T, N = [], [], []
    for y in ys:
        t2m, lat = gm(root, '2t', y)
        tsr, _ = gm(root, 'tsr', y)
        ttr, _ = gm(root, 'ttr', y)
        if t2m is None or tsr is None or ttr is None:
            continue
        w = np.cos(np.deg2rad(lat))
        Y.append(y)
        T.append(float(np.average(t2m.mean(axis=0).mean(axis=1), weights=w)) - 273.15)
        N.append(float(np.average(((tsr + ttr) / ACC).mean(axis=0).mean(axis=1), weights=w)))
    return np.array(Y), np.array(T), np.array(N)


def gregory(T, N):
    """slope [W/m2/K], the N=0 crossing temperature, and the slope's standard error.

    The crossing is returned only when the slope is significantly negative; otherwise
    it is unbounded and reporting it would invent an equilibrium.
    """
    if len(T) < 8:
        return np.nan, np.nan, np.nan
    b, a = np.polyfit(T, N, 1)
    res = N - (a + b * T)
    se = np.sqrt((res ** 2).sum() / (len(T) - 2) / ((T - T.mean()) ** 2).sum())
    teq = -a / b if (b < 0 and abs(b) > 1.96 * se) else np.nan
    return b, teq, se


# the late window is the same calendar years for every arm
LATE_START = min(max(sorted(int(f.split('_')[-1].split('-')[0])
                            for f in os.listdir(f'{r}/outdata/oifs')
                            if f.startswith('atm_remapped_1m_2t_') and f.endswith('.nc')))
                 for _l, r, _c in ARMS if r) - 14

fig, ax = plt.subplots(figsize=(9.4, 6.2))
ax.axhline(0, color='#2ca02c', lw=1.2, ls='--', zorder=1)
rows = []
for label, root, col in ARMS:
    if root is None:
        continue
    Y, T, N = series(root)
    k = Y > TRANSIENT_END
    Y, T, N = Y[k], T[k], N[k]
    if not len(Y):
        continue
    sc = ax.scatter(T, N, c=Y, cmap='viridis', s=26, edgecolor=col, linewidth=0.9,
                    zorder=3, label=f'{label}  {len(Y)} yr')
    half = Y >= LATE_START
    b_all, teq_all, se_all = gregory(T, N)
    b_late, teq_late, se_late = gregory(T[half], N[half])
    xs = np.linspace(T.min() - 0.05, T.max() + 0.35, 20)
    if np.isfinite(b_late):
        a_late = np.polyfit(T[half], N[half], 1)[1]
        ax.plot(xs, a_late + b_late * xs, color=col, lw=2.0, alpha=0.9, zorder=4)
    rows.append((label, len(Y), T[-8:].mean(), N[-8:].mean(),
                 b_all, se_all, teq_all, b_late, se_late, teq_late))

cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label('model year')
ax.set_xlabel('global mean T2m  [$^\\circ$C]')
ax.set_ylabel('global net TOA  [W m$^{-2}$]')
ax.set_title('Gregory view of the 1850 arms, transient years excluded\n'
             f'line = fit over {LATE_START}+, the same years for every arm',
             fontsize=11)
ax.legend(loc='lower left', fontsize=9)
ax.grid(alpha=0.25, lw=0.5)
fig.tight_layout()
p = f'{OUT}/gregory_1850_arms.png'
fig.savefig(p, dpi=160)
print(f'wrote {p}\n')

print(f'  late window: {LATE_START}+ , identical for every arm\n')
print(f'  {"arm":22s} {"yr":>3s} {"T last8":>8s} {"N last8":>8s} '
      f'{"slope all":>13s} {"slope late":>13s} {"T at N=0":>10s}')
for (label, n, tl, nl, ba, sa, ta, bl, sl, tel) in rows:
    f = lambda b, s: f'{b:+.3f}+-{s:.3f}' if np.isfinite(b) else '     -'
    te = f'{tel:8.2f} C' if np.isfinite(tel) else ' refused'
    print(f'  {label:22s} {n:3d} {tl:8.3f} {nl:+8.3f} {f(ba,sa):>13s} '
          f'{f(bl,sl):>13s} {te:>10s}')
print('\n  "refused" = the slope is not significantly negative, so the N=0 crossing is')
print('  unbounded.  It is left blank rather than extrapolated.')
