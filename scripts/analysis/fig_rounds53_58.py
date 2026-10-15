"""Report figures for rounds 53 to 58.

plots/arms_2169_ice.png       the one-change arms from the 2169 restart against the control and its second
                              realisation (ctl2): Antarctic February and September ice volume, Arctic September
                              volume, T2m 60-90N; 3-year running means, 2170 (cold-started atmosphere) left out
plots/spinup_v350cc_series.png  the CMIP7 spin-up 1850-1984: T2m, net TOA, sea-ice volume, iceberg melt and
                              the September area of deep convection south of 50S; 10-year running means

Usage: fig_rounds53_58.py   (module analysis-toolbox/python-04.2026)
"""
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
D = f'{REPO}/data/clim'

# ---- arms from the 2169 restart -------------------------------------------------------------------------
runs = [('tke_albsn082', 'control', 'k', 2.2, '-'), ('ctl2', 'control, second realisation', '0.55', 2.2, '-'),
        ('icb', 'icebergs', '#1f5fa8', 1.6, '-'), ('cdoi10', 'ice-ocean drag 0.010', '#c0392b', 1.6, '-'),
        ('wamz0', 'wave-model roughness', '#7fb069', 1.1, '--'), ('cpres15', 'ice strength c 15', '#b07aa1', 1.1, '--'),
        ('z0form', 'air-ice roughness (form)', '#e6a23c', 1.1, ':')]
panels = [('sivols_m02', 1e-3, '(a) Antarctic ice volume, February', '$10^3$ km$^3$'),
          ('sivols_m09', 1e-3, '(b) Antarctic ice volume, September', '$10^3$ km$^3$'),
          ('sivoln_m09', 1e-3, '(c) Arctic ice volume, September', '$10^3$ km$^3$'),
          ('t2m_6090n', 1.0, '(d) T2m 60-90N', '$^\\circ$C')]
fig, axs = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
for ax, (col, sc, title, unit) in zip(axs.ravel(), panels):
    for tag, lab, c, lw, ls in runs:
        s = pd.read_csv(f'{D}/series_{tag}.csv').set_index('year')[col] * sc
        s = s[(s.index >= 2171) & (s.index <= 2209)].rolling(3, center=True, min_periods=2).mean()
        ax.plot(s.index, s.values, color=c, lw=lw, ls=ls, label=lab)
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold')
    ax.set_ylabel(unit); ax.grid(alpha=0.3)
axs[1, 0].set_xlabel('model year'); axs[1, 1].set_xlabel('model year')
axs[0, 0].legend(fontsize=8.5, frameon=False, ncol=1, loc='upper right')
fig.suptitle('One-change arms from the 2169 restart (TKE only, hemispheric GM, albsn 0.82), 3-year running means',
             fontsize=12, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{REPO}/plots/arms_2169_ice.png', dpi=160)
plt.close(fig)

# ---- the CMIP7 spin-up ---------------------------------------------------------------------------------
s = pd.read_csv(f'{D}/series_spinup_v350cc.csv').set_index('year')
r = lambda x: x.rolling(10, center=True, min_periods=5).mean()
fig, axs = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
ax = axs[0, 0]
ax.plot(s.index, s['t2m_glob'], color='0.7', lw=0.8); ax.plot(s.index, r(s['t2m_glob']), color='k', lw=2)
ax.axhline(13.98, color='#c0392b', ls='--', lw=1.2, label='tuning line at the same age (1975-84)')
ax.set_title('(a) T2m global', loc='left', fontsize=11, fontweight='bold'); ax.set_ylabel('$^\\circ$C'); ax.legend(fontsize=8.5, frameon=False)
ax = axs[0, 1]
ax.plot(s.index, s['toa_net'], color='0.7', lw=0.8); ax.plot(s.index, r(s['toa_net']), color='k', lw=2)
ax.axhline(0, color='0.4', lw=0.8)
ax.set_title('(b) net TOA', loc='left', fontsize=11, fontweight='bold'); ax.set_ylabel('W m$^{-2}$')
ax = axs[1, 0]
for col, lab, c in (('siv_nh_sep', 'Arctic, September', '#1f5fa8'), ('siv_nh_mar', 'Arctic, March', '#9cc3e6'),
                    ('siv_sh_sep', 'Antarctic, September', '#c0392b'), ('siv_sh_feb', 'Antarctic, February x 5', '#e6a2a2')):
    x = s[col] * (5 if 'feb' in col else 1)
    ax.plot(s.index, r(x), color=c, lw=2, label=lab)
ax.set_title('(c) sea-ice volume', loc='left', fontsize=11, fontweight='bold'); ax.set_ylabel('$10^3$ km$^3$')
ax.legend(fontsize=8.5, frameon=False, ncol=2); ax.set_xlabel('model year')
ax = axs[1, 1]
ax.plot(s.index, r(s['icb_melt_Gt']), color='#1f5fa8', lw=2, label='iceberg melt (Gt yr$^{-1}$)')
ax.set_ylabel('Gt yr$^{-1}$', color='#1f5fa8')
ax2 = ax.twinx()
ax2.plot(s.index, r(s['so_mld1000_sep_Mkm2']), color='#c0392b', lw=2, label='September MLD > 1000 m south of 50S')
ax2.set_ylabel('$10^6$ km$^2$', color='#c0392b')
ax.axvline(1900, color='0.4', lw=0.8, ls=':')
ax.text(1901, ax.get_ylim()[0], ' iceberg grounding\n mode 1 and cap', fontsize=8, va='bottom', color='0.3')
ax.set_title('(d) iceberg melt and Southern Ocean deep convection', loc='left', fontsize=11, fontweight='bold')
ax.set_xlabel('model year')
for a in axs.ravel():
    a.grid(alpha=0.3)
fig.suptitle('AWI-ESM3-CC-LR-CMIP7-Spinup, 1850-1984 (thin: annual, thick: 10-year running mean)', fontsize=12, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{REPO}/plots/spinup_v350cc_series.png', dpi=160)
plt.close(fig)
print('written: plots/arms_2169_ice.png plots/spinup_v350cc_series.png')
