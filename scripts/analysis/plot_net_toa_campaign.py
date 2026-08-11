"""The radiation objective, brought up to date: coupled through round 11, and AMIP.

Regenerates campaign_net_toa_by_year_with10.png with the 11 series added, and produces
the AMIP companion the coupled figure has always implied but never shown.

WHY TWO FIGURES.  The coupled plot carries a horizontal line labelled "atmosphere-only
floor (AMIP, +0.67)" -- a single number standing in for the entire atmosphere-only
campaign.  That was fine when the AMIP runs all sat near it.  They no longer do: the
round-27 L-series moves global net TOA from +0.64 to between +1.08 and -0.88, i.e. it
straddles the goal and the floor is not a floor any more.  So the second panel plots the
AMIP arms as time series on the same axis convention, and the "floor" line in the coupled
figure is relabelled as what it actually is -- the AMIP CONTROL, one arm among many.

DESIGN.  Generations are distinguished by WEIGHT, not just hue: 09 and 10 are drawn thin
and muted because they are settled history, the 11 series bold because it is the live
question.  End-of-line direct labels with a minimum separation, no legend box -- with this
many series a legend forces a colour lookup for every line, and the lines end at different
heights anyway.  Colour follows the run, never its rank.

  net TOA = (tsr + ttr) / accumulation      both positive DOWN, ttr negative

VALIDATION.  The 09 values come from a cached CSV while everything else is recomputed from
raw output, so 09C is recomputed and compared against its cache before anything is
plotted.  IFS TOA fluxes are accumulated J/m2 over the output step; getting that wrong
once made an M-series pass come out 3600x too large.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import csv, numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

BASE = ('/work/ab0246/a270092/postprocessing/'
        'investigation_awiesm3_high_lat_cold_bias_global_toa_positive')
OUT, DAT = f'{BASE}/plots/', f'{BASE}/data/'
STEP, GOAL = 3600.0, 0.0
SURF, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983'
RT270 = '/work/bb1469/a270270/runtime/awiesm3-v3.4'
RT092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
RTAMIP = '/work/bb1469/a270092/runtime/oifsamip-cy48'

CACHED = {'080a (CRUNCEP veg init)': ('080a  baseline, old sea ice', '#c9c3bb', 0),
          '09A = baseline + newSeaIce': ('09A  baseline + newSeaIce', '#bdd8c9', 0),
          '09B = 06T + newSeaIce': ('09B  06T + newSeaIce', '#e3d3a8', 0),
          '09C = 06V + newSeaIce': ('09C  06V + newSeaIce', '#e8bccd', 0)}
# (label, root, colour, emphasis)  emphasis 0 = settled history, 1 = live
COUPLED = [
    ('10A  G4', f'{RT270}/Tuning_test_10A_06V_G4_CRUNCEP_plus_CERES_init_newSeaIce', '#b7a6d4', 0),
    ('10B  G4 + tanh depletion', f'{RT270}/Tuning_test_10B_06V_G4_snowDepletion_CRUNCEP_plus_CERES_init_newSeaIce', '#a8c6dd', 0),
    ('110  09C + IFSsoilT', f'{RT270}/Tuning_test_110Baseline_09C_useIFSsoiltemp_CRUNCEPinit_newSeaIce', '#8a8983', 1),
    ('11A  + G4 + D2b', f'{RT270}/Tuning_test_11A_06V_G4_useIFSsoiltemp_CRUNCEP_plus_CERES_init_newSeaIce', '#1baf7a', 1),
    ('11B  + tanh snow (falsified)', f'{RT270}/Tuning_test_11B_06V_G4_snowDepletion_useIFSsoiltemp_CRUNCEP_plus_CERES_init_newSeaIce', '#eb6834', 1),
    ('11D  + fitted snow, SWEMIN 30', f'{RT270}/Tuning_test_11D_G4_fitted_snow_depletion_useIFSsoiltemp_CRUNCEP_plus_CERES_init_newSeaIce', '#eda100', 1),
    ('11E  + SWEMIN 15 + K1', f'{RT092}/Tuning_test_11E_swemin15_K1', '#7b4fb5', 1),
    ('11F  + DMS 166 (cancelled)', f'{RT092}/Tuning_test_11F_dmsccn166', '#1f7ab5', 1),
]
AMIP = [
    ('control  (defaults + RVICE .16)', 'amip_pi_base', '#8a8983', 1),
    ('P5  adopted stack', 'amip_P5_swemin15', '#1baf7a', 1),
    ('S4  + INPPMIN 50k', 'amip_S4_inppmin50000', '#eda100', 1),
    ('LX1  + RSNOWLIN2', 'amip_LX1_long', '#7b4fb5', 1),
    ('LX3  + DMS + INPPMIN', 'amip_LX3_long', '#eb6834', 1),
    ('LY2  + overlap 0.10', 'amip_LY2_long', '#1f7ab5', 1),
]
VALIDATE = ('09C = 06V + newSeaIce', f'{RT092}/Tuning_test_09C_06V_CRUNCEPinit_newSeaIce')

print(__doc__)


def net_toa(root, y0=1300, y1=2000):
    out, d = {}, f'{root}/outdata/oifs'
    if not os.path.isdir(d):
        return out
    for f in sorted(os.listdir(d)):
        if not f.startswith('atm_remapped_1m_tsr_') or not f.endswith('.nc'):
            continue
        try:
            yr = int(f.split('_')[-1].split('-')[0])
        except ValueError:
            continue
        if not (y0 <= yr <= y1):
            continue
        g = f.replace('_tsr_', '_ttr_')
        if not os.path.exists(f'{d}/{g}'):
            continue
        with xr.open_dataset(f'{d}/{f}', decode_times=False) as a:
            tsr, lat = a['tsr'].values, a['lat'].values
        with xr.open_dataset(f'{d}/{g}', decode_times=False) as b:
            ttr = b['ttr'].values
        net = (tsr + ttr) / STEP
        ww = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], net.shape[1:])
        out[yr] = float(np.average(net.mean(axis=0), weights=ww))
    return out


# ---------------------------------------------------------------- validation
cached = defaultdict(list)
csvp = f'{DAT}campaign_net_toa_by_year.csv'
if os.path.exists(csvp):
    for r in csv.DictReader(open(csvp)):
        try:
            cached[r['label']].append((int(r['year']), float(r['netTOA'])))
        except ValueError:
            pass
lab, root = VALIDATE
chk, ref = net_toa(root), dict(cached.get(lab, []))
common = sorted(set(chk) & set(ref))
if common:
    dif = np.array([chk[y] - ref[y] for y in common])
    ok = np.abs(dif).max() <= 0.02
    print(f'VALIDATION {lab}: {len(common)} yr, max |diff| {np.abs(dif).max():.4f} '
          f'-> {"agrees" if ok else "*** DISAGREES, units/timestep suspect ***"}')
else:
    print('VALIDATION: no overlap; proceeding UNVALIDATED')

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.linewidth': 0.6,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.labelcolor': INK2, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF})


def sm(y, k=5):
    y = np.asarray(y); r = np.full(len(y), np.nan)
    for i in range(len(y)):
        r[i] = y[max(0, i - k // 2):min(len(y), i + k // 2 + 1)].mean()
    return r


def draw(series, title, ref_line, ref_lab, out_png, xlim, ylim, mins=0.17):
    fig, ax = plt.subplots(figsize=(11.8, 5.2))
    fig.subplots_adjust(left=0.075, right=0.70, top=0.915, bottom=0.11)
    ax.axhline(GOAL, color=INK2, lw=1.1, ls=(0, (5, 3)), zorder=2)
    if ref_line is not None:
        ax.axhline(ref_line, color='#2a78d6', lw=1.3, ls=(0, (1.5, 2)), zorder=2)
        ax.axhspan(GOAL, ref_line, color='#2a78d6', alpha=0.06, lw=0, zorder=1)
        ax.text(xlim[0] + 0.6, ref_line + 0.05, ref_lab, fontsize=8.2, color='#2a78d6',
                zorder=10, bbox=dict(facecolor=SURF, edgecolor='none', alpha=.78, pad=1.5))
    ax.text(xlim[0] + 0.6, GOAL + 0.05, 'piControl goal (0.0)', fontsize=8.2, color=INK2,
            zorder=10, bbox=dict(facecolor=SURF, edgecolor='none', alpha=.78, pad=1.5))
    ends = []
    for disp, col, x, y, emph in series:
        lw_raw, lw_sm, al = (0.9, 2.3, .30) if emph else (0.6, 1.3, .18)
        ax.plot(x, y, color=col, lw=lw_raw, alpha=al, zorder=3)
        ys = sm(y)
        ax.plot(x, ys, color=col, lw=lw_sm, solid_capstyle='round',
                alpha=1.0 if emph else 0.75, zorder=4 + emph)
        ax.plot(x[-1], ys[-1], 'o', ms=5.5 if emph else 3.8, color=col, mec=SURF,
                mew=1.2, zorder=5 + emph)
        ends.append([ys[-1], disp, col, emph])
    ends.sort(key=lambda e: e[0])
    for i in range(1, len(ends)):
        if ends[i][0] - ends[i - 1][0] < mins:
            ends[i][0] = ends[i - 1][0] + mins
    for yv, disp, col, emph in ends:
        ax.annotate(disp, xy=(1.02, yv), xycoords=('axes fraction', 'data'),
                    fontsize=8.4 if emph else 7.8, va='center', color=col,
                    alpha=1.0 if emph else 0.75, annotation_clip=False,
                    fontweight='bold' if emph else 'normal')
    ax.set_xlabel('model year')
    ax.set_ylabel('global net TOA imbalance  [W m$^{-2}$]')
    ax.grid(axis='y', color=MUTED, alpha=0.22, lw=0.5); ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    fig.suptitle(title, x=0.008, ha='left', fontsize=11.5, fontweight='bold',
                 color=INK, y=0.995)
    fig.savefig(out_png, dpi=170)
    print('  saved:', out_png)


# ---------------------------------------------------------------- coupled figure
series = []
for key, (disp, col, emph) in CACHED.items():
    xy = sorted(cached.get(key, []))
    if xy:
        series.append((disp, col, np.array([a for a, _ in xy]),
                       np.array([b for _, b in xy]), emph))
print('\nCOUPLED:')
for disp, root, col, emph in COUPLED:
    d = net_toa(root)
    if not d:
        print(f'  {disp}: no output'); continue
    xs = np.array(sorted(d)); ys = np.array([d[y] for y in xs])
    series.append((disp, col, xs, ys, emph))
    print(f'  {disp:32s} {xs[0]}-{xs[-1]}  last decade {ys[-10:].mean():+.3f}')
draw(series, 'The radiation objective — coupled, rounds 09 to 11',
     0.6684, 'AMIP control (+0.67)', OUT + 'campaign_net_toa_by_year_with11.png',
     (1349.5, 1400.5), (-1.2, 2.5))

# ---------------------------------------------------------------- AMIP figure
print('\nAMIP:')
aser = []
for disp, run, col, emph in AMIP:
    d = net_toa(f'{RTAMIP}/{run}', 1860, 1930)
    if not d:
        print(f'  {disp}: no output'); continue
    xs = np.array(sorted(d)); ys = np.array([d[y] for y in xs])
    aser.append((disp, col, xs, ys, emph))
    ev = [d[y] for y in xs if 1872 <= y <= 1915]
    print(f'  {disp:34s} {xs[0]}-{xs[-1]}  1872-1915 mean {np.mean(ev):+.3f}')
draw(aser, 'The radiation objective — atmosphere-only (AMIP), the round-27 candidates',
     None, None, OUT + 'campaign_net_toa_by_year_amip.png',
     (1869.5, 1918.5), (-1.6, 2.2), mins=0.14)
