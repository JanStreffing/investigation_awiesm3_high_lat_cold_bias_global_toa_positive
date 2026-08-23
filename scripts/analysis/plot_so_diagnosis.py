"""Figure: the Southern Ocean has two separate problems, in two different seasons.

(a) The SH ice seasonal cycle is too large: too much ice at the September maximum, far
    too little through February-May, and the March regrowth is absent entirely.
(b) The T2m warm anomaly over the ice zone peaks in MAM, the REFREEZE season, not in the
    melt season -- which no shortwave mechanism explains.
(c) A +0.9 K warm anomaly at exactly 100 m, directly under a neutral surface layer and
    above COLD water at 500-1000 m.  That is the depth the autumn mixed layer reaches.

Numbers from so_seaice_and_albedo.py, coupled_sb_eval.py and so_subsurface_vs_phc3.py;
this only draws them, so the four must be kept in step.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

MON = ['J','F','M','A','M','J','J','A','S','O','N','D']
obs_a = [3.62, 2.23, 3.10, 5.82, 8.97, 11.86, 14.19, 15.70, 16.27, 15.63, 12.91, np.nan]
m11p  = [4.00, 1.87, 1.85, 4.10, 7.81, 11.38, 14.24, 16.30, 17.57, 16.67, 13.52, 8.62]
m11r  = [3.19, 1.31, 1.30, 3.37, 6.95, 10.66, 13.74, 16.07, 17.05, 16.18, 13.30, 8.34]

fig, ax = plt.subplots(1, 3, figsize=(13.2, 3.9))

a = ax[0]
x = np.arange(12)
a.plot(x, obs_a, 'k-o', lw=2, ms=4, label='OSI-SAF', zorder=3)
a.plot(x, m11p, '-s', color='#2166ac', lw=1.6, ms=3.5, label='11P')
a.plot(x, m11r, '-^', color='#66bd63', lw=1.6, ms=3.5, label='11R')
a.fill_between(x[1:5], np.array(obs_a[1:5]), np.array(m11r[1:5]),
               color='#d73027', alpha=0.18, zorder=1)
a.annotate('Feb–May deficit\n(no March regrowth)', xy=(2.6, 2.4), xytext=(3.2, 8.3),
           fontsize=8, color='#b2182b',
           arrowprops=dict(arrowstyle='->', color='#b2182b', lw=1.0))
a.set_xticks(x); a.set_xticklabels(MON, fontsize=8)
a.set_ylabel('SH sea ice area [$10^6$ km$^2$]')
a.set_title('(a)  seasonal cycle too large', fontsize=10, loc='left')
a.legend(fontsize=8, frameon=False, loc='upper left'); a.grid(alpha=0.25)

b = ax[1]
seas = ['DJF', 'MAM', 'JJA', 'SON']
w = 0.36
b11p = [-0.40, 1.21, 1.03, 0.48]
b11r = [-0.21, 1.69, 1.53, 0.66]
xs = np.arange(4)
b.bar(xs - w/2, b11p, w, color='#2166ac', label='11P')
b.bar(xs + w/2, b11r, w, color='#66bd63', label='11R')
b.axhline(0, color='0.4', lw=0.8)
b.annotate('peaks in the\nREFREEZE season', xy=(1.0, 1.69), xytext=(1.25, 1.15),
           fontsize=8, color='#b2182b',
           arrowprops=dict(arrowstyle='->', color='#b2182b', lw=1.0))
b.set_xticks(xs); b.set_xticklabels(seas, fontsize=9)
b.set_ylabel('T$_{2m}$ bias, 60–90$^\\circ$S ocean [K]')
b.set_title('(b)  warm anomaly is an autumn problem', fontsize=10, loc='left')
b.legend(fontsize=8, frameon=False, loc='upper left'); b.grid(alpha=0.25, axis='y')

c = ax[2]
dep = np.array([0, 50, 100, 200, 500, 1000, 2000])
d11p = np.array([-0.25, -0.15, +0.87, +0.22, -0.40, -0.42, -0.12])
d11r = np.array([-0.17, -0.07, +0.93, +0.29, -0.33, -0.34, -0.09])
c.plot(d11p, dep, '-s', color='#2166ac', lw=1.6, ms=4, label='11P')
c.plot(d11r, dep, '-^', color='#66bd63', lw=1.6, ms=4, label='11R')
c.axvline(0, color='0.4', lw=0.8)
c.axhspan(60, 160, color='#d73027', alpha=0.13, zorder=0)
c.annotate('+0.9 K at 100 m:\nthe autumn mixed\nlayer reaches this', xy=(0.93, 100),
           xytext=(0.05, 430), fontsize=8, color='#b2182b',
           arrowprops=dict(arrowstyle='->', color='#b2182b', lw=1.0))
c.invert_yaxis(); c.set_yscale('symlog', linthresh=200)
c.set_yticks([0, 50, 100, 200, 500, 1000, 2000])
c.set_yticklabels(['0', '50', '100', '200', '500', '1000', '2000'], fontsize=8)
c.set_xlabel('model $-$ PHC3 [K]'); c.set_ylabel('depth [m]')
c.set_title('(c)  90–60S subsurface', fontsize=10, loc='left')
c.legend(fontsize=8, frameon=False, loc='lower right'); c.grid(alpha=0.25)

fig.tight_layout()
out = 'report/run_figures/so_two_problems.png'
fig.savefig(out, dpi=170)
print('wrote', out)
