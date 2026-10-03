"""Final evaluation of the four vertical-mixing arms, 2020-2039.

All four branch from PICAL_ccnice at 2019-12-31 and differ from it and from each other in
ONE namelist line, mix_scheme, so this is a paired comparison on identical years from an
identical restart - which this campaign has rarely had.

  PICAL_ccnice        KPP                     (control, FESOM's own KPP)
  PICAL_cvKPP         cvmix_KPP               (the null: same scheme, other implementation)
  PICAL_cvTKE         cvmix_TKE
  PICAL_cvTKEIDEMIX   cvmix_TKE+cvmix_IDEMIX

Two figures:
  plots/mixing_eval_volume.png  sea-ice volume, the four seasonal extremes, per year
  plots/mixing_eval_state.png   the state metrics, as distance from the observational target

OBSERVATIONAL REFERENCES, and the epoch caveat that governs both figures: GIOMAS and
OSI-SAF are 1989-1999 (pre-2000, the better pre-industrial analogue) and CERES is present
day, while these runs are pre-industrial.  In the NH the offset is large and the model
should sit ABOVE the reference; in the SH Antarctic ice carried no significant trend over
that record, so PD ~ PI is defensible and the panel reads near face value.

Usage:  python3 scripts/figures/mixing_scheme_eval.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings, textwrap
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

REPO = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive'
R    = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
Y0, Y1 = 2020, 2039
SURF, INK, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#d8d7d2'

# validated categorical order, slots 1-4
ARMS = [('PICAL_ccnice',      'KPP  (control)',            '#2a78d6'),
        ('PICAL_cvKPP',       'cvmix_KPP',                 '#eda100'),
        ('PICAL_cvTKE',       'cvmix_TKE',                 '#1baf7a'),
        ('PICAL_cvTKEIDEMIX', 'cvmix_TKE + cvmix_IDEMIX',  '#eb6834')]

def vol(arm, var, y):
    p = f'{R}/{arm}/outdata/fesom/{var}.fesom.{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p, decode_times=False) as d:
        a = np.squeeze(np.asarray(d[var].values, float))/1e3     # units 1e9 m3 -> 10^3 km3
    return a if a.size == 12 else None

# ---------------------------------------------------------------- figure 1: volume
PANELS = [('NH April  (maximum)', 'sivoln', 3, 31.19),
          ('NH September  (minimum)', 'sivoln', 8, 11.96),
          ('SH September  (maximum)', 'sivols', 8, 19.89),
          ('SH February  (minimum)', 'sivols', 1, 2.07)]

fig, axes = plt.subplots(2, 2, figsize=(12.4, 7.4))
fig.patch.set_facecolor(SURF)
for ax, (title, var, mon, obs) in zip(axes.ravel(), PANELS):
    ax.set_facecolor(SURF)
    ax.axhline(obs, color=MUTED, lw=1.4, ls='--', zorder=1)
    ax.annotate(f'observed {obs:.1f}', (Y0+0.2, obs), textcoords='offset points',
                xytext=(0, 4), fontsize=8, color=MUTED)
    for arm, lab, col in ARMS:
        ys, vs = [], []
        for y in range(Y0, Y1+1):
            a = vol(arm, var, y)
            if a is not None: ys.append(y); vs.append(a[mon])
        if not ys: continue
        last = (arm == 'PICAL_cvTKEIDEMIX')
        ax.plot(ys, vs, color=col, lw=2.4 if last else 1.9, marker='o', ms=3.4,
                mec=SURF, mew=0.7, label=lab, zorder=4 if last else 3)
    ax.set_title(title, fontsize=10.5, color=INK, loc='left', fontweight='bold')
    ax.grid(alpha=0.18, lw=0.6); ax.tick_params(labelsize=8, colors=MUTED)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'): ax.spines[s].set_color(GRID)
axes[0, 0].set_ylabel('volume  [10$^3$ km$^3$]', fontsize=9, color=MUTED)
axes[1, 0].set_ylabel('volume  [10$^3$ km$^3$]', fontsize=9, color=MUTED)
axes[0, 0].legend(frameon=False, fontsize=9, loc='upper right', labelcolor=MUTED, ncol=2)
fig.suptitle('Vertical mixing schemes: sea-ice volume, 2020-2039', fontsize=13,
             color=INK, x=0.008, ha='left', y=0.985)
note = ('All four branch from PICAL_ccnice at 2019-12-31 and differ in one namelist line. '
        'Observed = GIOMAS 1989-1999. The runs are pre-industrial and the reference is not: '
        'in the NH the model should sit ABOVE the dashed line, in the SH Antarctic ice carried '
        'no significant trend over that record so the panel reads near face value.')
fig.text(0.008, 0.005, '\n'.join(textwrap.wrap(note, 150)), fontsize=7.6, color=MUTED, va='bottom')
fig.tight_layout(rect=[0, 0.055, 1, 0.945])
p1 = f'{REPO}/plots/mixing_eval_volume.png'
fig.savefig(p1, dpi=160, facecolor=SURF); plt.close(fig); print('wrote', p1)

# ---------------------------------------------------------------- figure 2: state scorecard
# Each metric as distance from its observational target, signed so that NEGATIVE is always
# "too little / too cold" and positive "too much / too warm".  Plotted as distance so the
# arms are comparable across quantities with different units.
METRICS = [
    # label,                        target, values per arm in ARMS order, units
    ('SH Feb extent',                 3.00, [1.847, 1.659, 1.952, 2.263], '10$^6$ km$^2$'),
    ('SH Sep extent',                18.50, [21.696, 21.736, 20.921, 21.313], '10$^6$ km$^2$'),
    ('SH Sep volume',                19.89, [10.99, 11.75, 14.23, 14.22], '10$^3$ km$^3$'),
    ('T2m 60-90S',                   -0.75, [0.757, 0.824, 0.664, 0.569], 'K'),
    ('SW CRE 45-65S',               -68.02, [-64.463, -64.839, -65.760, -65.325], 'W m$^{-2}$'),
    ('NH Apr volume',                31.19, [33.64, 36.70, 35.17, 38.04], '10$^3$ km$^3$'),
    ('NH Mar extent',                15.00, [18.070, 19.161, 17.780, 18.446], '10$^6$ km$^2$'),
    ('T2m global',                   -0.85, [-0.461, -0.491, -0.317, -0.505], 'K'),
    ('net TOA',                       0.00, [0.215, -0.044, -0.436, -0.141], 'W m$^{-2}$'),
]

fig, ax = plt.subplots(figsize=(11.2, 6.4)); fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
n = len(METRICS); h = 0.19
for j, (arm, lab, col) in enumerate(ARMS):
    xs = [(vals[j]-tgt)/abs(tgt) if tgt != 0 else (vals[j]-tgt)
          for _, tgt, vals, _ in METRICS]
    ys = [n-1-i + (j-1.5)*h for i in range(n)]
    ax.barh(ys, xs, height=h*0.88, color=col, label=lab, zorder=3)
ax.axvline(0, color=INK, lw=1.4, zorder=4)
ax.set_yticks(range(n))
ax.set_yticklabels([f'{m[0]}\n{MUTED and ""}target {m[1]:g} {m[3]}' for m in METRICS][::-1],
                   fontsize=8.5, color=MUTED)
ax.set_xlabel('fractional distance from the observational target   (0 = on target)',
              fontsize=9, color=MUTED)
ax.tick_params(axis='x', labelsize=8, colors=MUTED)
ax.grid(axis='x', alpha=0.18, lw=0.6)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
for s in ('left', 'bottom'): ax.spines[s].set_color(GRID)
ax.legend(frameon=False, fontsize=9, loc='lower right', labelcolor=MUTED)
fig.suptitle('Vertical mixing schemes: distance from the observational target, 2020-2039',
             fontsize=13, color=INK, x=0.008, ha='left', y=0.985)
note2 = ('Bars are (model - target)/|target|, so shorter is better and the sign says which way. '
         'net TOA is an absolute difference, its target being zero. The pre-industrial epoch '
         'offset is NOT removed: the NH rows and T2m global are expected to overshoot their '
         'satellite-era references, the SH rows much less so.')
fig.text(0.008, 0.005, '\n'.join(textwrap.wrap(note2, 140)), fontsize=7.6, color=MUTED, va='bottom')
fig.tight_layout(rect=[0, 0.075, 1, 0.945])
p2 = f'{REPO}/plots/mixing_eval_state.png'
fig.savefig(p2, dpi=160, facecolor=SURF); plt.close(fig); print('wrote', p2)
