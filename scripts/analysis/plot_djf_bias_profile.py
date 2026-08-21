"""Figure: the DJF land cold bias is surface-confined, and the skin is not the problem.

Panel (a) shows the model-minus-ERA5 DJF bias over 60-90N land against height for the
three coupled arms.  Panel (b) splits the near-surface error into the air-side inversion
bias (T925-T2m) and the skin-side one (T2m-Tskt).  Together they are the argument that
retired the surface-side levers: the bias decays to nothing by 850 hPa, so it is a
redistribution problem, but the skin is already right, so the error is in the layer
between the screen and 925 hPa.

Numbers come from djf_bias_vertical_structure.py; this only draws them, so the two must
be kept in step.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 60-90N land, DJF, model minus ERA5.  From djf_bias_vertical_structure.py, 2026-08-21.
ARMS = {
    '11G inppmin50k': dict(prof=[-6.345, -4.373, -2.815, -2.123], lowinv=3.520, sfc=0.346),
    '11N LX4 1850':   dict(prof=[-4.491, -2.095, -1.329, -1.085], lowinv=2.884, sfc=0.276),
    '11P LX4 1990':   dict(prof=[-3.627, -2.100, -0.876, -0.418], lowinv=2.809, sfc=0.110),
}
Y = [1013, 1000, 925, 850]          # 2 m plotted just below 1000 hPa
LAB = ['T$_{2m}$', '1000', '925', '850']
COL = {'11G inppmin50k': '#b2182b', '11N LX4 1850': '#ef8a62', '11P LX4 1990': '#2166ac'}

fig, ax = plt.subplots(1, 2, figsize=(9.4, 4.0), gridspec_kw=dict(width_ratios=[1.25, 1]))

a = ax[0]
for k, v in ARMS.items():
    a.plot(v['prof'], Y, 'o-', color=COL[k], lw=1.8, ms=5, label=k)
a.axvline(0, color='0.4', lw=0.8)
a.invert_yaxis()
a.set_yticks(Y); a.set_yticklabels(LAB)
a.set_xlabel('model $-$ ERA5  [K]')
a.set_ylabel('level')
a.set_title('(a)  DJF bias by height, 60–90$^\\circ$N land', fontsize=10, loc='left')
a.legend(fontsize=8, frameon=False, loc='lower left')
a.grid(alpha=0.25)
a.annotate('bias decays to\nnothing aloft', xy=(-0.42, 850), xytext=(-3.4, 883),
           fontsize=8, color='#2166ac',
           arrowprops=dict(arrowstyle='->', color='#2166ac', lw=1.0))

b = ax[1]
x = np.arange(len(ARMS)); w = 0.36
lo = [v['lowinv'] for v in ARMS.values()]
sf = [v['sfc'] for v in ARMS.values()]
b.bar(x - w/2, lo, w, color='#2166ac', label='air side:  T$_{925}-$T$_{2m}$')
b.bar(x + w/2, sf, w, color='#bbbbbb', label='skin side: T$_{2m}-$T$_{skt}$')
for xi, v in zip(x - w/2, lo):
    b.text(xi, v + 0.08, f'{v:+.2f}', ha='center', fontsize=8, color='#2166ac')
for xi, v in zip(x + w/2, sf):
    b.text(xi, v + 0.08, f'{v:+.2f}', ha='center', fontsize=8, color='0.35')
b.axhline(0, color='0.4', lw=0.8)
b.set_xticks(x); b.set_xticklabels([k.split()[0] for k in ARMS], fontsize=9)
b.set_ylabel('inversion bias  [K]')
b.set_ylim(-0.2, 4.1)
b.set_title('(b)  where the near-surface error sits', fontsize=10, loc='left')
b.legend(fontsize=8, frameon=False, loc='upper right')
b.grid(alpha=0.25, axis='y')

fig.tight_layout()
out = 'report/run_figures/djf_bias_profile.png'
fig.savefig(out, dpi=170)
print('wrote', out)
