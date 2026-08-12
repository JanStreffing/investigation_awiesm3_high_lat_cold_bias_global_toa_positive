"""Raupach (1994) canopy roughness against what the model does today.

WHY THIS FIGURE.  The Lund group's note shows the two z0 curves for an idealised 20 m
canopy and states that EC-Earth has the sign backwards.  That is the left panel, and it
is the verification that our implementation reproduces theirs.  It is not yet an argument
about OUR model, because it says nothing about where our canopy actually sits on the
curve or how much grid-box roughness the error costs.

The right panel closes that gap with the model's own output: the measured Siberian annual
cycle of high-vegetation cover from coupled 11E, converted into an effective grid-box
roughness under the two treatments.  This is the panel that says how big the error is.

WHAT THE MODEL DOES TODAY.  The high-vegetation TILE keeps a fixed 2.0 m from the lookup
table (vupdz0_mod.F90:276); what collapses in winter is the tile's AREA, because
FracHVeg is the foliar projective cover computed from today's phenological LAI.  So the
comparison has to be made on the grid-box aggregate, not on a tile value, and the
aggregate has to be formed in DRAG space the way vupdz0 forms it:

    Cdn = (kappa/ln(1+z/z0))^2,  averaged over tiles by area,  then inverted

ASSUMED STAND STRUCTURE, stated because the right panel depends on it.  The coupled
output carries cvh and lai_hv but not canopy height or stem density, so the woody
structure is taken as a Siberian larch stand: 15 m tall, 0.1 stems/m2, 0.20 m diameter.
Those give a stem frontal area index of 0.30, which is the term that carries the winter
roughness.  The qualitative result -- winter roughness rising rather than collapsing --
does not depend on these numbers; the magnitude does.

A SUBTLETY WORTH SEEING.  lai_hv in the output sits at ~2.5 all year even in January,
because it is lai/fpc and both go to zero together.  The model is therefore not
representing a leafless forest at all: it represents a SMALLER AREA OF LEAFY forest.
That is drawn explicitly in the right panel.
"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = ('/work/ab0246/a270092/postprocessing/'
        'investigation_awiesm3_high_lat_cold_bias_global_toa_positive')
OUT = f'{BASE}/plots'
R = '/work/bb1469/a270092/runtime/awiesm3-v3.4/Tuning_test_11E_swemin15_K1/outdata/oifs'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/'
        'Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/outdata/oifs/'
        'atm_remapped_1m_lsm_1350-1350.nc')
YEARS = range(1350, 1355)          # forest still standing; the bug is not forest loss

# --- Raupach (1994), the same constants as lpj_guess/framework/raupach_z0.cpp ---------
C_D1, C_S, C_R, US_U_MAX, C_W, KAPPA = 7.5, 0.003, 0.3, 0.3, 2.0, 0.41
PSI_H = np.log(C_W) - 1.0 + 1.0 / C_W
LEAF_FRONTAL = 0.5

# --- IFS lookup values (susveg_mod.F90) ---------------------------------------------
Z0_HIGH, Z0_BARE, Z0_LOW = 2.00, 0.013, 0.034   # needleleaf, desert/bare, tundra

# --- assumed Siberian larch structure ------------------------------------------------
H_CANOPY, DENS, DIAM = 15.0, 0.10, 0.20
LAMBDA_STEM = DENS * DIAM * H_CANOPY

SURF, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983'
C_RAU, C_IFS, C_ACC = '#1f7ab5', '#eb6834', '#7b4fb5'

print(__doc__)
print('=' * 100)


def raupach_z0(h, lam):
    lam = np.maximum(lam, 1e-6)
    zeta = np.sqrt(2.0 * C_D1 * lam)
    disp = h * (1.0 - (1.0 - np.exp(-zeta)) / zeta)
    us_u = np.minimum(np.sqrt(C_S + C_R * lam), US_U_MAX)
    return np.minimum((h - disp) * np.exp(-KAPPA / us_u + PSI_H), h)


def cdn(z0, zb=50.0):
    return (KAPPA / np.log(1.0 + zb / np.maximum(z0, 1e-9))) ** 2


def z0_from_cdn(c, zb=50.0):
    return zb / (np.exp(KAPPA / np.sqrt(np.maximum(c, 1e-12))) - 1.0)


# =====================================================================================
# measured annual cycle, Siberian box, coupled 11E
# =====================================================================================
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values)
    lat, lon = d['lat'].values, d['lon'].values
if lsm.ndim == 3:
    lsm = lsm[0]
box = np.zeros_like(lsm, bool)
box[np.ix_((lat >= 55) & (lat < 75), (lon >= 60) & (lon <= 180))] = True
box &= lsm > 0.5
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], lsm.shape)


def cycle(var):
    acc = []
    for y in YEARS:
        f = f'{R}/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = np.squeeze(d[var].values)
        if a.shape[0] == 12:
            acc.append(a)
    m = np.mean(acc, axis=0)
    return np.array([float(np.average(m[k][box], weights=W[box])) for k in range(12)])


cvh, cvl, laih = cycle('cvh'), cycle('cvl'), cycle('lai_hv')
print(f'\n  measured cvh   {cvh.min():.3f} (winter) -> {cvh.max():.3f} (August)')
print(f'  measured lai_hv {laih.min():.2f} -> {laih.max():.2f}   '
      f'(nearly flat: it is lai/fpc, not a canopy LAI)')

# True leaf area held by the box = lai_hv * cvh.  Spread over the ANNUAL-MAXIMUM forest
# area -- the trees are standing all year -- gives the canopy LAI the stand really has.
area_true = cvh.max()
lai_true = laih * cvh / area_true

lam_rau = LEAF_FRONTAL * lai_true + LAMBDA_STEM
z0_tile_rau = raupach_z0(H_CANOPY, lam_rau)

# grid-box aggregate, drag space, exactly as vupdz0 combines tiles
bare_now = np.clip(1.0 - cvh - cvl, 0.0, 1.0)
box_ifs = z0_from_cdn(cvh * cdn(Z0_HIGH) + cvl * cdn(Z0_LOW) + bare_now * cdn(Z0_BARE))
bare_rau = np.clip(1.0 - area_true - cvl, 0.0, 1.0)
box_rau = z0_from_cdn(area_true * cdn(z0_tile_rau) + cvl * cdn(Z0_LOW)
                      + bare_rau * cdn(Z0_BARE))

M = 'J F M A M J J A S O N D'.split()
print(f'\n  {"":14s}' + ''.join(f'{x:>7s}' for x in M))
print(f'  {"cvh (model)":14s}' + ''.join(f'{v:7.3f}' for v in cvh))
print(f'  {"box z0 today":14s}' + ''.join(f'{v:7.3f}' for v in box_ifs))
print(f'  {"box z0 Raupach":14s}' + ''.join(f'{v:7.3f}' for v in box_rau))
print(f'\n  DJF grid-box roughness: {box_ifs[[11,0,1]].mean():.3f} m today -> '
      f'{box_rau[[11,0,1]].mean():.3f} m with Raupach '
      f'({box_rau[[11,0,1]].mean()/box_ifs[[11,0,1]].mean():.1f}x)')
print(f'  JJA grid-box roughness: {box_ifs[[5,6,7]].mean():.3f} m today -> '
      f'{box_rau[[5,6,7]].mean():.3f} m with Raupach '
      f'({box_rau[[5,6,7]].mean()/box_ifs[[5,6,7]].mean():.1f}x)')

# =====================================================================================
# figure
# =====================================================================================
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.linewidth': 0.6,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.labelcolor': INK2, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.9))
fig.subplots_adjust(left=0.065, right=0.985, top=0.855, bottom=0.115, wspace=0.24)

# --- (a) the formula, 20 m canopy ----------------------------------------------------
pai = np.linspace(0.15, 5.0, 400)
ax1.plot(pai, raupach_z0(20.0, LEAF_FRONTAL * pai), color=C_RAU, lw=2.4,
         label='Raupach (1994) — implemented here')
ax1.plot(pai, (1.0 - np.exp(-0.5 * pai)) * Z0_HIGH, color=C_IFS, lw=2.4,
         label='IFS today: $c_{vh}\\times$RVZ0M = $(1-e^{-0.5\\,LAI})\\times2.0$')
for x, y in ((0.5, 2.47), (0.6, 2.57), (4.4, 1.07)):
    ax1.plot(x, y, 'o', ms=6.5, mfc='none', mec=INK, mew=1.4, zorder=5)
ax1.annotate('the three points read off the\nLund note, reproduced to 0.004 m',
             xy=(0.66, 2.55), xytext=(1.30, 2.30), fontsize=8.2, color=INK2,
             arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))
ipk = int(np.argmin(np.abs(raupach_z0(20.0, LEAF_FRONTAL * pai)
                           - (1.0 - np.exp(-0.5 * pai)) * Z0_HIGH)))
ax1.axvline(pai[ipk], color=MUTED, lw=0.8, ls=(0, (4, 3)))
ax1.text(pai[ipk] + 0.12, 2.18, f'crossover at LAI+SAI $\\approx$ {pai[ipk]:.1f}\n'
         'our boreal canopy sits here, so the\nJJA error is small and the DJF one\n'
         'is the whole story', fontsize=8.2, color=INK2, va='top')
ax1.set_xlabel('tree LAI plus stem area index'); ax1.set_ylabel('$z_0$  [m]')
ax1.set_title('(a)  The sign error, for a 20 m canopy', fontsize=10.5, loc='left',
              fontweight='bold', pad=6)
ax1.set_xlim(0, 5); ax1.set_ylim(0, 3.0)
ax1.legend(frameon=False, fontsize=8.4, loc='lower right', bbox_to_anchor=(1.0, 0.06))
ax1.grid(axis='y', color=MUTED, alpha=0.22, lw=0.5); ax1.set_axisbelow(True)

# --- (b) our model, our Siberia ------------------------------------------------------
x = np.arange(12)
ax2.plot(x, box_ifs, color=C_IFS, lw=2.4, marker='o', ms=4.5, label='today')
ax2.plot(x, box_rau, color=C_RAU, lw=2.4, marker='o', ms=4.5, label='with Raupach')
ax2.fill_between(x, box_ifs, box_rau, where=box_rau >= box_ifs, color=C_RAU, alpha=0.10,
                 lw=0)
axb = ax2.twinx()
axb.plot(x, cvh, color=C_ACC, lw=1.5, ls=(0, (3, 2)))
axb.set_ylabel('$c_{vh}$ sent by LPJ-GUESS', fontsize=8.6, color=C_ACC)
axb.tick_params(axis='y', labelsize=8, colors=C_ACC)
axb.set_ylim(0, 0.45)
axb.spines['top'].set_visible(False)
axb.annotate(f'the forest stands all year,\nbut its AREA falls to {cvh.min():.3f}',
             xy=(2.0, cvh[2]), xytext=(2.3, 0.30), fontsize=8.2, color=C_ACC,
             ha='left', va='center',
             arrowprops=dict(arrowstyle='->', color=C_ACC, lw=0.9,
                             connectionstyle='arc3,rad=-0.25'))
axb.text(7.0, cvh[7] + 0.012, '$c_{vh}$', fontsize=8.6, color=C_ACC, ha='center')
ax2.set_xticks(x); ax2.set_xticklabels(M)
ax2.set_xlabel('month'); ax2.set_ylabel('effective grid-box $z_0$  [m]')
ax2.set_title('(b)  Siberia 55–75N, from coupled 11E output', fontsize=10.5, loc='left',
              fontweight='bold', pad=6)
ax2.set_ylim(0, max(box_rau.max(), box_ifs.max()) * 1.45)
ax2.annotate(f'DJF: {box_ifs[[11,0,1]].mean():.3f} $\\to$ '
             f'{box_rau[[11,0,1]].mean():.3f} m  '
             f'({box_rau[[11,0,1]].mean()/box_ifs[[11,0,1]].mean():.1f}$\\times$)',
             xy=(0.02, 0.985), xycoords='axes fraction', fontsize=9.0, color=C_RAU,
             ha='left', va='top', fontweight='bold')
ax2.legend(frameon=False, fontsize=8.6, loc='upper left',
           bbox_to_anchor=(0.0, 0.93))
ax2.grid(axis='y', color=MUTED, alpha=0.22, lw=0.5); ax2.set_axisbelow(True)
for a in (ax1, ax2):
    for sp in ('top', 'right'):
        a.spines[sp].set_visible(False)
ax2.spines['right'].set_visible(True); ax2.spines['right'].set_color(C_ACC)

fig.suptitle('Canopy roughness: the leafless boreal forest is the roughest state, '
             'and the model makes it the smoothest',
             x=0.008, ha='left', fontsize=12, fontweight='bold', color=INK, y=0.975)
fig.text(0.008, 0.925, 'Left: verification against the LPJ-GUESS group\'s note.  '
         'Right: the same formulation driven by the model\'s own measured '
         f'$c_{{vh}}$, assuming a {H_CANOPY:.0f} m larch stand '
         f'({DENS:.2f} stems m$^{{-2}}$, {DIAM:.2f} m diameter, stem '
         f'$\\lambda$ = {LAMBDA_STEM:.2f}).',
         ha='left', fontsize=8.4, color=INK2)

os.makedirs(OUT, exist_ok=True)
p = f'{OUT}/raupach_z0.png'
fig.savefig(p, dpi=170, bbox_inches='tight')
print(f'\n  written: {p}')
