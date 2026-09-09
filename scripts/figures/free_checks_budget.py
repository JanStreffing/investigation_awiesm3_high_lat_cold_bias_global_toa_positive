"""Report figure: the drift is linear, and the bands that accumulate are fed laterally."""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
d = np.load(os.path.join(REPO, 'data', 'free_checks_11X.npz'), allow_pickle=True)
th, so, yrs, bud, bands = d['theta'], d['soga'], d['years'], d['budget'], d['bands']

fig, ax = plt.subplots(1, 2, figsize=(15, 5.2))

# (a) linearity of the drift
a = ax[0]
a.plot(yrs, th, 'o-', color='#c0392b', ms=4, lw=1.4, label=r'$\theta_{\rm global}$')
c = np.polyfit(yrs[5:], th[5:], 1)   # fit the ramp, not the initial adjustment
a.plot(yrs[5:], np.polyval(c, yrs[5:]), 'k--', lw=1.4,
       label=f'linear fit 1355-89, {c[0]*10:+.4f} K/decade')
for i in range(0, 40, 10):
    a.axvspan(yrs[i], yrs[i] + 9, color='0.9' if (i // 10) % 2 else 'white', zorder=0)
a.set_xlabel('year')
a.set_ylabel(r'global mean potential temperature [$^\circ$C]')
a.set_title('(a) the drift does not decelerate: +1.461, +1.433, +1.472 W m$^{-2}$\n'
            'across the three decade steps', fontsize=11, loc='left')
a.grid(alpha=0.3)
a.legend(fontsize=9)
a2 = a.twinx()
a2.plot(yrs, so, '-', color='#2980b9', lw=1.2)
a2.set_ylabel('global mean salinity [psu]', color='#2980b9')
a2.tick_params(axis='y', labelcolor='#2980b9')
a2.set_ylim(so.mean() - 0.002, so.mean() + 0.002)
a2.text(0.97, 0.06, 'soga: $-2.5\\times10^{-5}$ psu in 30 yr (salt is conserved)',
        fontsize=8, color='#2980b9', ha='right', transform=a2.transAxes)
a.annotate('flat until ~1355,\nthen straight', xy=(1353, th[3]), xytext=(1357, 3.700),
           fontsize=8, color='0.35',
           arrowprops=dict(arrowstyle='->', color='0.5', lw=0.8))

# (b) band budget
b = ax[1]
y = np.arange(len(bands))
b.barh(y - 0.22, bud[:, 0], height=0.2, color='#e08214', label='into ocean at the surface')
b.barh(y, bud[:, 1], height=0.2, color='#4a90c0', label='lateral convergence')
b.barh(y + 0.22, bud[:, 2], height=0.2, color='#7b2d26', label='kept (heat content gain)')
b.axvline(0, color='0.3', lw=1.0)
b.set_yticks(y)
b.set_yticklabels(bands, fontsize=9)
b.invert_yaxis()
b.set_xlabel('W m$^{-2}$ of band area')
b.set_title('(b) only the tropics and the subantarctic gain heat at the surface;\n'
            'every band that accumulates is fed laterally', fontsize=11, loc='left')
b.grid(alpha=0.3, axis='x')
b.legend(fontsize=9, loc='lower right')

fig.tight_layout()
out = os.path.join(REPO, 'report', 'run_figures', 'free_checks_budget_11X.png')
fig.savefig(out, dpi=140, bbox_inches='tight')
print('wrote', out)
