"""Two figures: what Raupach roughness does, and why it cannot work here.

FIG 1 (mechanism).  z0 as a function of the roughness density lambda, evaluated with the
implementation's OWN constants, for two canopy heights: the 2.89 m scrub the coupled model
actually grows in Siberia, and an 18 m boreal larch stand.  The operating points and the
IFS lookup value are marked.  The message the eye should get without reading a word: the
model sits on the steep RISING branch far below the peak, and below the IFS line it
replaces -- so switching Raupach on makes the winter surface smoother.

FIG 2 (measurement).  The 11J-11I seasonal momentum response with 95 % CI, three metrics
in small multiples (native units -- never a dual axis), DJF and JJA side by side, with the
pre-registered DJF sign marked.  The message: summer moved and is significant; winter did
not move at all, and what movement there is has the wrong sign.

Palette: dataviz reference categorical slots 1 (blue) and 2 (orange), validated
  node scripts/validate_palette.js "#2a78d6,#eb6834" --mode light  -> ALL CHECKS PASS
  worst adjacent CVD dE 24.7 (protan), normal-vision 33.6, contrast >= 3:1.
Identity is never colour-alone: both series carry direct labels and distinct markers.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '..', '..', 'report', 'plots')
os.makedirs(OUT, exist_ok=True)

BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, INK2, MUTED = '#1a1a19', '#4a4a46', '#8a8a82'
GRID = '#e3e3df'

# ---- the implementation's own constants, framework/raupach_z0.cpp -------------
C_D1, C_S, C_R, US_U_MAX, KAPPA = 7.5, 0.003, 0.3, 0.3, 0.41
PSI_H = 0.193                      # ln(c_w) - 1 + 1/c_w, c_w = 2.0
MIN_H, MIN_LAM = 0.5, 1.0e-4


def raupach_z0(h, lam):
    if h < MIN_H or lam < MIN_LAM:
        return np.nan
    z = np.sqrt(2.0 * C_D1 * lam)
    d = h * (1.0 - (1.0 - np.exp(-z)) / z)
    usu = min(np.sqrt(C_S + C_R * lam), US_U_MAX)
    return min((h - d) * np.exp(-KAPPA / usu + PSI_H), h)


# ---- measured in 11J's Siberian box, 728 cells, 1370-79 ----------------------
H_MODEL, LAM_MODEL = 2.89, 0.0092          # 2 cm stems, 0.161 indiv/m2
H_REAL, LAM_REAL = 18.0, 0.12              # a real larch stand, leaf-off
IFS_DJF, IFS_JJA, BARE = 0.058, 0.190, 0.013

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED,
                     'axes.labelcolor': INK, 'text.color': INK,
                     'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.linewidth': 0.8, 'figure.dpi': 200})

# =============================== FIGURE 1 =====================================
fig, ax = plt.subplots(figsize=(6.9, 4.3))
lam = np.logspace(-3.1, 0.3, 400)

for h, col, lab, ls in ((H_REAL, BLUE, f'real boreal stand, $h={H_REAL:.0f}$ m', '-'),
                        (H_MODEL, ORANGE, f'coupled model, $h={H_MODEL:.2f}$ m', '-')):
    z = np.array([raupach_z0(h, l) for l in lam])
    ax.plot(lam, z, ls, color=col, lw=2.0, zorder=3)
    # direct label, so identity is never colour-alone
    i = np.argmin(np.abs(lam - 0.55))
    ax.annotate(lab, xy=(lam[i], z[i]), xytext=(6, -2), textcoords='offset points',
                color=col, fontsize=8.5, fontweight='bold', va='center')

# IFS lookup values it replaces
for v, txt in ((IFS_DJF, 'IFS lookup, DJF  0.058 m'), (IFS_JJA, 'IFS lookup, JJA  0.190 m')):
    ax.axhline(v, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.text(1.1e-3, v * 1.13, txt, color=INK2, fontsize=7.6, va='bottom')

# operating points, with a 2 px surface ring so they read over the lines
for (h, l, col, lab, dx, dy) in (
        (H_MODEL, LAM_MODEL, ORANGE,
         f'11J winter\n$\\lambda_{{stem}}={LAM_MODEL:.4f}$\n$z_0={raupach_z0(H_MODEL, LAM_MODEL):.4f}$ m', 14, -34),
        (H_REAL, LAM_REAL, BLUE,
         f'real stand, leaf-off\n$\\lambda_{{stem}}={LAM_REAL:.2f}$\n$z_0={raupach_z0(H_REAL, LAM_REAL):.2f}$ m', -96, 6)):
    z = raupach_z0(h, l)
    ax.plot([l], [z], 'o', ms=9, color=col, mec='white', mew=2.0, zorder=5)
    ax.annotate(lab, xy=(l, z), xytext=(dx, dy), textcoords='offset points',
                fontsize=7.8, color=INK,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=col, lw=0.9, alpha=0.95),
                arrowprops=dict(arrowstyle='-', color=col, lw=0.9), zorder=6)

# the finding, stated on the plot
ax.annotate('', xy=(LAM_MODEL, raupach_z0(H_MODEL, LAM_MODEL)), xytext=(LAM_MODEL, IFS_DJF),
            arrowprops=dict(arrowstyle='-|>', color=ORANGE, lw=1.8, shrinkA=1, shrinkB=4),
            zorder=4)
ax.text(LAM_MODEL * 0.62, np.sqrt(IFS_DJF * raupach_z0(H_MODEL, LAM_MODEL)),
        r'$\mathbf{4.4\times}$' + '\nsmoother', color=ORANGE, fontsize=8.4,
        fontweight='bold', va='center', ha='right', linespacing=1.15)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e-3, 2.0); ax.set_ylim(5e-3, 6)
ax.set_xlabel(r'roughness density  $\lambda$  (frontal area index)  [-]')
ax.set_ylabel(r'roughness length  $z_0$  [m]')
ax.set_title('Raupach (1994) $z_0$: the winter mechanism needs a canopy the model has not got',
             fontsize=9.8, fontweight='bold', pad=9, loc='left')
ax.grid(True, which='major', color=GRID, lw=0.6, zorder=0)
ax.grid(True, which='minor', color=GRID, lw=0.3, alpha=0.6, zorder=0)
ax.set_axisbelow(True)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
fig.tight_layout()
f1 = os.path.join(OUT, 'raupach_z0_curve.png')
fig.savefig(f1, bbox_inches='tight')
plt.close(fig)
print('wrote', f1)

# =============================== FIGURE 2 =====================================
# 11J - 11I, Siberian land, 30 paired years 1350-79 (raupach_column_response.py)
M = [('$|\\tau|$  [N m$^{-2}$]',
      dict(DJF=(-0.0039, -0.0087, 0.0010, 0.113), JJA=(-0.0047, -0.0085, -0.0009, 0.016))),
     ('10 m wind  [m s$^{-1}$]',
      dict(DJF=(+0.0629, -0.0473, 0.1731, 0.253), JJA=(+0.0816, 0.0054, 0.1578, 0.037))),
     ('skin $-$ 2 m air  [K]',
      dict(DJF=(-0.0042, -0.0353, 0.0270, 0.787), JJA=(+0.1041, 0.0615, 0.1468, 0.00003)))]
SEASON_Y = {'DJF': 1.0, 'JJA': 0.0}
SEASON_C = {'DJF': BLUE, 'JJA': ORANGE}

fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.5))
for ax, (name, vals) in zip(axes, M):
    ax.axvline(0, color=MUTED, lw=1.0, zorder=1)
    for s, (d, lo, hi, p) in vals.items():
        y = SEASON_Y[s]
        c = SEASON_C[s]
        ax.plot([lo, hi], [y, y], '-', color=c, lw=2.0, solid_capstyle='round', zorder=3)
        ax.plot([d], [y], 'o', ms=8, color=c, mec='white', mew=1.8, zorder=4)
        sig = p < 0.05
        # p below the bar, so the predicted-sign arrow row above stays clear
        ax.text(hi, y - 0.30, ('$p=%.5f$' % p) if p < 1e-3 else ('$p=%.3f$' % p),
                fontsize=7.0, color=(INK if sig else MUTED), ha='right', va='top',
                fontweight=('bold' if sig else 'normal'))
    ax.set_yticks([SEASON_Y['DJF'], SEASON_Y['JJA']])
    ax.set_yticklabels(['DJF', 'JJA'], fontweight='bold')
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.75, 1.55)
    ax.set_title(name, fontsize=8.6, pad=5)
    ax.grid(True, axis='x', color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis='x', labelsize=7.4)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))

# The pre-registered winter prediction: |tau| UP, wind DOWN.  Drawn in AXES
# coordinates so it can never fall outside the data range, and only on the two
# panels where the mechanical sign is unambiguous (skin-air is a consequence, not
# a prediction of the roughness formula).
for ax_i, frac_from, frac_to in ((0, 0.62, 0.92), (1, 0.38, 0.08)):
    a = axes[ax_i]
    a.annotate('', xy=(frac_to, 0.90), xytext=(frac_from, 0.90),
               xycoords='axes fraction', textcoords='axes fraction',
               arrowprops=dict(arrowstyle='-|>', color=MUTED, lw=1.2))
    a.text((frac_from + frac_to) / 2, 0.955, 'predicted DJF', transform=a.transAxes,
           fontsize=6.8, color=MUTED, ha='center', va='bottom')

fig.suptitle('11J $-$ 11I, Siberian land, 30 paired years: summer responds, winter does not',
             fontsize=9.8, fontweight='bold', x=0.005, ha='left', y=1.06)
fig.text(0.005, -0.10,
         'Points are paired mean differences, bars 95 % CI. Winter was predicted rougher '
         '($|\\tau|$ up, wind down) and instead has the same signs as\nsummer. Only '
         'JJA skin$-$air survives Bonferroni over all 14 tests.',
         fontsize=7.2, color=INK2, ha='left', linespacing=1.35)
fig.tight_layout()
f2 = os.path.join(OUT, 'raupach_momentum_response.png')
fig.savefig(f2, bbox_inches='tight')
plt.close(fig)
print('wrote', f2)
