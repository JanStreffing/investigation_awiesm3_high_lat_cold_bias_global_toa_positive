"""Figure: the stable-BL trade is linear in dose, and coupling changes its terms.

Panel (a) plots what each arm PAYS (tropical SW CRE) against what it BUYS (Arctic DJF
inversion removed).  SB1 and SB2 sit on one line through the origin -- 0.88 and 0.87
W/m2 per K -- so halving the RSBLB dose halves both, and the trade cannot be dodged by
tuning gentler.  SB3, the mixing-length route, sits far off that line at 3.43.

Panel (b) is why the coupled pair was run.  In AMIP the inversion weakens 40 % by warming
the screen and 60 % by cooling the air above, and the screen response is not significant.
Coupled, with sea ice free to respond, it is essentially all screen.

Numbers from sb_series_eval.py and coupled_sb_eval.py; this only draws them.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# (benefit = -Arctic DJF inversion delta, cost = tropical SW CRE delta), 1874-1901 AMIP
AMIP = {'SB1 $b$=3': (0.318, 0.281, '#4393c3'),
        'SB2 $b$=2': (0.624, 0.541, '#2166ac'),
        'SB3 $\\lambda_{min}$=120': (0.614, 2.106, '#b2182b')}
CPL = {'11Q (1850)': (0.680, 0.486, '#1a9850'),
       '11R (1990)': (0.809, 0.540, '#66bd63')}

fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.2))

a = ax[0]
xs = np.linspace(0, 0.95, 20)
a.plot(xs, 0.875 * xs, '--', color='0.55', lw=1.2, zorder=1)
a.text(0.60, 0.60, 'RSBLB frontier\n0.87 W m$^{-2}$ per K', fontsize=8, color='0.4')
# hand-placed offsets: SB2, 11Q and 11R sit close together on the frontier
OFF = {'SB1 $b$=3': (8, -14), 'SB2 $b$=2': (-16, -20),
       'SB3 $\\lambda_{min}$=120': (8, 8), '11Q (1850)': (10, 6), '11R (1990)': (10, -4)}
for k, (b_, c_, col) in {**AMIP, **CPL}.items():
    mk = 'o' if k.startswith('SB') else 's'
    a.scatter(b_, c_, s=95, color=col, marker=mk, zorder=3, edgecolor='w', lw=1.0)
    a.annotate(k, (b_, c_), textcoords='offset points', xytext=OFF[k],
               fontsize=8.5, color=col)
a.set_xlabel('benefit: Arctic DJF inversion removed [K]')
a.set_ylabel('cost: tropical SW CRE [W m$^{-2}$]')
a.set_title('(a)  the trade is linear, and SB3 is off it', fontsize=10, loc='left')
a.set_xlim(0, 1.0); a.set_ylim(0, 2.35)
a.grid(alpha=0.25)
a.scatter([], [], marker='o', color='0.4', label='AMIP')
a.scatter([], [], marker='s', color='0.4', label='coupled')
a.legend(fontsize=8, frameon=False, loc='upper left')

b = ax[1]
lab = ['AMIP SB2', 'coupled 11Q', 'coupled 11R']
screen = [0.285, 0.614, 0.674]          # DJF T2m 60-90N
aloft = [-0.428, -0.114, -0.020]        # DJF T925 60-90N
sig = [False, True, True]
x = np.arange(3); w = 0.36
b.bar(x - w/2, screen, w, color=['#bbbbbb', '#1a9850', '#66bd63'],
      label='screen  $\\Delta$T$_{2m}$')
b.bar(x + w/2, aloft, w, color='#2166ac', alpha=0.85, label='aloft  $\\Delta$T$_{925}$')
for xi, v, s in zip(x - w/2, screen, sig):
    b.text(xi, v + 0.03, f'{v:+.2f}' + ('*' if s else ' (ns)'), ha='center', fontsize=8)
for xi, v in zip(x + w/2, aloft):
    b.text(xi, v - 0.09, f'{v:+.2f}', ha='center', fontsize=8, color='#2166ac')
b.axhline(0, color='0.4', lw=0.8)
b.set_xticks(x); b.set_xticklabels(lab, fontsize=9)
b.set_ylabel('DJF change, 60–90$^\\circ$N land [K]')
b.set_ylim(-0.62, 0.95)
b.set_title('(b)  coupling moves the response to the screen', fontsize=10, loc='left')
b.legend(fontsize=8, frameon=False, loc='upper left')
b.grid(alpha=0.25, axis='y')

fig.tight_layout()
out = 'report/run_figures/sb_exchange_rate.png'
fig.savefig(out, dpi=170)
print('wrote', out)
