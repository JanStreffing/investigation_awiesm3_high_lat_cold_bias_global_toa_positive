"""The radiation objective, extended with the 10 series.

Extends the round-09 net-TOA figure (plot2.py, recovered 2026-07-28) with the coupled
10 series, which branches from the same point and runs 20 years longer (1350-1399
against 1350-1379).

WHAT THE 10 SERIES ADDS.  10A is 09C's configuration plus G4 (the boreal RVRSMIN
lever); 10B is 10A plus the tanh snow depletion -- the scheme that has since been
falsified outright (surfece.F90, and the AMIP N/O series).  So the pair shows what G4
bought in the coupled model, and what the broken depletion scheme cost it.

VALIDATION BUILT IN.  The 09-series values are read from the cached CSV, but the 10
series has to be computed from raw output, so the two could disagree through a units
or timestep error.  The script therefore RECOMPUTES 09C from its own output and
compares against the cached number before plotting anything.  IFS TOA fluxes are
accumulated J/m2 over the output step -- the same trap that made the first M-series
pass come out 3600x too large -- so this check is not optional.

    net TOA = (tsr + ttr) / step        both positive DOWN, ttr negative
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import csv, numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

BASE = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive'
OUT, DAT = f'{BASE}/report/plots/', f'{BASE}/data/'
STEP = 3600.0
AMIP_FLOOR = 0.6684
# TARGET CHANGED 2026-08-07: 0.0, not the inherited HR piControl -0.16.
# -0.16 was never derived -- it is whatever the HR piControl settled at, adopted as a
# matching target, and report.tex already flagged that it should be revisited.  A
# target away from zero is only justified if the atmosphere fails to conserve energy,
# because equilibrium requires SFC = TOA + S for a spurious atmospheric source S, so
# zero OCEAN drift would need TOA = -S.  Measured on 10A (1390-1399, global, using the
# release tool's own formula ssr+str+sshf+slhf - sf*3.3355e8, all /accumulation):
#     net TOA +0.854, net surface +0.877  ->  residual -0.023 W/m2
# The atmosphere conserves to 0.02, so "ocean not drifting" and "system not drifting"
# are the same condition and both put the target at 0.  Snow enthalpy is 0.966 W/m2
# and does nearly all the reconciliation; omitting it makes the surface look like
# +1.84 against a TOA of +0.85, which is the artefact that motivated -0.16-style fudges.
GOAL = 0.0

SURF, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983'
RT270 = '/work/bb1469/a270270/runtime/awiesm3-v3.4'
RT092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'

# cached 09-series labels -> display
CACHED = {'080a (CRUNCEP veg init)': ('080a  baseline, old sea ice', '#eb6834'),
          '09A = baseline + newSeaIce': ('09A  baseline + newSeaIce', '#1baf7a'),
          '09B = 06T + newSeaIce': ('09B  06T + newSeaIce', '#eda100'),
          '09C = 06V + newSeaIce': ('09C  06V + newSeaIce', '#e87ba4')}
# runs computed here
COMPUTE = [('10A  G4', f'{RT270}/Tuning_test_10A_06V_G4_CRUNCEP_plus_CERES_init_newSeaIce', '#7b4fb5'),
           ('10B  G4 + tanh depletion', f'{RT270}/Tuning_test_10B_06V_G4_snowDepletion_CRUNCEP_plus_CERES_init_newSeaIce', '#1f7ab5')]
VALIDATE = ('09C = 06V + newSeaIce', f'{RT092}/Tuning_test_09C_06V_CRUNCEPinit_newSeaIce')

print(__doc__)


def net_toa(root, y0=1300, y1=1500):
    """Annual global-mean net TOA [W/m2] per year."""
    out = {}
    d = f'{root}/outdata/oifs'
    if not os.path.isdir(d):
        return out
    for f in sorted(os.listdir(d)):
        if not f.startswith('atm_remapped_1m_tsr_') or not f.endswith('.nc'):
            continue
        yr = int(f.split('_')[-1].split('-')[0])
        if not (y0 <= yr <= y1):
            continue
        g = f.replace('_tsr_', '_ttr_')
        if not os.path.exists(f'{d}/{g}'):
            continue
        with xr.open_dataset(f'{d}/{f}', decode_times=False) as a:
            tsr = a['tsr'].values; lat = a['lat'].values
        with xr.open_dataset(f'{d}/{g}', decode_times=False) as b:
            ttr = b['ttr'].values
        net = (tsr + ttr) / STEP
        w = np.cos(np.deg2rad(lat))
        ww = np.broadcast_to(w[:, None], net.shape[1:])
        out[yr] = float(np.average(net.mean(axis=0), weights=ww))
    return out


# ---- validation ---------------------------------------------------------------
cached = defaultdict(list)
for r in csv.DictReader(open(f'{DAT}campaign_net_toa_by_year.csv')):
    try:
        cached[r['label']].append((int(r['year']), float(r['netTOA'])))
    except ValueError:
        pass

lab, root = VALIDATE
chk = net_toa(root)
ref = dict(cached[lab])
common = sorted(set(chk) & set(ref))
if common:
    dif = np.array([chk[y] - ref[y] for y in common])
    print(f'VALIDATION on {lab}: {len(common)} common years, '
          f'mean diff {dif.mean():+.4f}, max |diff| {np.abs(dif).max():.4f} W/m2')
    if np.abs(dif).max() > 0.02:
        print('  *** RECOMPUTATION DISAGREES WITH THE CACHE -- units/timestep suspect ***')
    else:
        print('  -> agrees; the same recipe is safe for the 10 series')
else:
    print(f'VALIDATION: no overlapping years for {lab}; proceeding UNVALIDATED')

# ---- the 10 series --------------------------------------------------------------
series = []
for name, _, col in [(k, None, v[1]) for k, v in CACHED.items()]:
    pass
for key, (disp, col) in CACHED.items():
    xy = sorted(cached[key])
    if xy:
        series.append((disp, col, np.array([a for a, _ in xy]), np.array([b for _, b in xy])))
for disp, root, col in COMPUTE:
    d = net_toa(root)
    if not d:
        print(f'  {disp}: no output'); continue
    xs = np.array(sorted(d)); ys = np.array([d[y] for y in xs])
    series.append((disp, col, xs, ys))
    print(f'  {disp}: {xs[0]}-{xs[-1]}, last decade {ys[-10:].mean():+.3f} W/m2')

# ---- plot -----------------------------------------------------------------------
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.linewidth': 0.6,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.labelcolor': INK2, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF})


def sm(y, k=5):
    y = np.asarray(y); r = np.full(len(y), np.nan)
    for i in range(len(y)):
        r[i] = y[max(0, i - k // 2):min(len(y), i + k // 2 + 1)].mean()
    return r


fig, ax = plt.subplots(figsize=(11.4, 5.0))
# top margin allows TWO subtitle lines without colliding with the title
fig.subplots_adjust(left=0.08, right=0.715, top=0.92, bottom=0.11)
ax.axhline(GOAL, color=INK2, lw=1.1, ls=(0, (5, 3)), zorder=2)
ax.axhline(AMIP_FLOOR, color='#2a78d6', lw=1.4, ls=(0, (1.5, 2)), zorder=2)
ax.axhspan(GOAL, AMIP_FLOOR, color='#2a78d6', alpha=0.06, lw=0, zorder=1)

ends = []
for disp, col, x, y in series:
    ax.plot(x, y, color=col, lw=0.9, alpha=0.28, zorder=3)
    ys = sm(y)
    ax.plot(x, ys, color=col, lw=2.2, solid_capstyle='round', zorder=4)
    ax.plot(x[-1], ys[-1], 'o', ms=5.5, color=col, mec=SURF, mew=1.3, zorder=5)
    ends.append([ys[-1], disp, col])
ends.sort(key=lambda e: e[0])
MINS = 0.17
for i in range(1, len(ends)):
    if ends[i][0] - ends[i - 1][0] < MINS:
        ends[i][0] = ends[i - 1][0] + MINS
for yv, disp, col in ends:
    ax.annotate(disp, xy=(1.02, yv), xycoords=('axes fraction', 'data'), fontsize=8.4,
                va='center', color=col, annotation_clip=False)

ax.text(1351, AMIP_FLOOR + 0.06, 'atmosphere-only floor (AMIP, +0.67)', fontsize=8.2,
        color='#2a78d6', zorder=10,
        bbox=dict(facecolor=SURF, edgecolor='none', alpha=0.78, pad=1.5))
ax.text(1351, GOAL + 0.06, 'piControl goal (0.0)', fontsize=8.2, color=INK2, zorder=10,
        bbox=dict(facecolor=SURF, edgecolor='none', alpha=0.78, pad=1.5))
ax.set_xlabel('model year'); ax.set_ylabel('global net TOA imbalance  [W m$^{-2}$]')
ax.grid(axis='y', color=MUTED, alpha=0.22, lw=0.5); ax.set_axisbelow(True)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
ax.set_xlim(1349.5, 1400.5); ax.set_ylim(-1.0, 2.5)
fig.suptitle('The radiation objective — rounds 09 and 10 vs the atmosphere-only floor',
             x=0.008, ha='left', fontsize=11.5, fontweight='bold', color=INK, y=0.995)
# subtitle removed on request -- the goal-change rationale lives in the header
# comment above, in RUNS_AND_PARAMETERS.md and in report.tex, not on the figure.
out = OUT + 'campaign_net_toa_by_year_with10.png'
fig.savefig(out, dpi=170)
print('\nSaved:', out)

with open(DAT + 'campaign_net_toa_by_year_with10.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['series', 'year', 'netTOA'])
    for disp, _, x, y in series:
        for a, b in zip(x, y):
            w.writerow([disp, int(a), round(float(b), 4)])
print('Saved:', DAT + 'campaign_net_toa_by_year_with10.csv')
