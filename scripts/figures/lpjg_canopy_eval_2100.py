"""Figures for the 2100 LPJ-GUESS merge (canopy work) evaluation on albedo.

1. plots/lpjg_canopy_siberia_veg.png   Siberian box vegetation: larch survival from the
   CRUNCEP state (crunveg vs 080a, 11G) and the production line across the 2100 switch.
2. plots/lpjg_canopy_climate_ice.png   T2m and NH sea-ice volume on the production line
   across the switch, with crunveg (same restarts, forested land) beside it.

Inputs are the CSVs written by boreal_veg_trajectory.py and climate_annual_series.py.
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

V, C = 'data/veg_traj', 'data/clim'
COL = {'ccnice': '#2a78d6', 'crunveg': '#eb6834', '080a': '#1baf7a', '11G': '#eda100'}
LS = {'ccnice': '-', 'crunveg': '-', '080a': '--', '11G': ':'}
INK, MUTED, GRID = '#1a1a19', '#6b6a64', '#e4e3dd'
plt.rcParams.update({'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.6, 'axes.spines.top': False,
                     'axes.spines.right': False, 'font.size': 9, 'lines.linewidth': 2})


def veg():
    ccn = pd.concat([pd.read_csv(f'{V}/ccnice_lev_sib.csv'), pd.read_csv(f'{V}/ccnice_alb_sib.csv')])
    crun = pd.read_csv(f'{V}/crunveg_sib.csv')
    old = {k: pd.read_csv(f'{V}/{k}_sib.csv') for k in ('080a', '11G')}
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))

    a = ax[0]
    for k, d in [('crunveg', crun)] + list(old.items()):
        t = d.year - d.year.min()
        a.plot(t, d.fpc_BNS, c=COL[k], ls=LS[k],
               label={'crunveg': 'PICAL_crunveg (merged LPJ-GUESS, 2100+)',
                      '080a': '080a (pre-merge, 1350+)', '11G': '11G (pre-merge, 1350+)'}[k])
    a.set(xlabel='years since the CRUNCEP land state', ylabel='BNS (larch) FPC',
          title='a  Larch from the CRUNCEP start', ylim=(0, None))
    a.legend(fontsize=8, frameon=False)

    a = ax[1]
    for k, d in [('crunveg', crun)] + list(old.items()):
        t = d.year - d.year.min()
        a.plot(t, d.dens_BNS, c=COL[k], ls=LS[k], label=k)
    a.set(xlabel='years since the CRUNCEP land state', ylabel='BNS density (ind m$^{-2}$)',
          title='b  Larch density: mortality outruns establishment in 080a/11G', ylim=(0, None))
    a.legend(fontsize=8, frameon=False)

    a = ax[2]
    s = ccn[ccn.year >= 2040]
    a.plot(s.year, s.fpc_TREEFPC, c=COL['ccnice'], label='PICAL_ccnice tree FPC')
    a.plot(s.year, s.fpc_GRASSFPC, c=COL['ccnice'], ls='--', label='PICAL_ccnice grass FPC')
    a.plot(crun.year, crun.fpc_TREEFPC, c=COL['crunveg'], label='PICAL_crunveg tree FPC')
    a.plot(crun.year, crun.fpc_GRASSFPC, c=COL['crunveg'], ls='--', label='PICAL_crunveg grass FPC')
    a.axvline(2100, c=MUTED, lw=1)
    a.text(2101, 0.34, 'LPJ-GUESS merge,\nmove to albedo', color=MUTED, fontsize=8, va='top')
    a.set(xlabel='model year', ylabel='FPC', title='c  Production line: no forest returns',
          ylim=(0, 0.36))
    a.legend(fontsize=7.5, frameon=False, loc='lower left')
    fig.suptitle('Siberian box 55-75N / 60-180E, cos(lat)-weighted LPJ-GUESS annual output',
                 color=INK)
    fig.tight_layout()
    fig.savefig('plots/lpjg_canopy_siberia_veg.png', dpi=140)


def climate():
    ccn = pd.concat([pd.read_csv(f'{C}/ccnice_lev.csv'), pd.read_csv(f'{C}/ccnice_alb.csv')])
    crun = pd.read_csv(f'{C}/crunveg.csv')
    ccn = ccn[ccn.year >= 2040]
    panels = [('t2m_glob', 'global T2m (°C)', None),
              ('t2m_6090n', 'T2m 60-90N (°C)', None),
              ('sivoln_m04', 'NH April ice volume (10$^3$ km$^3$)', 31.19),
              ('sivoln_m09', 'NH September ice volume (10$^3$ km$^3$)', 11.96)]
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.9))
    for a, (k, lab, obs) in zip(ax, panels):
        sc = 1e-3 if k.startswith('si') else 1
        a.plot(ccn.year, ccn[k] * sc, c=COL['ccnice'], label='PICAL_ccnice (grassland)')
        a.plot(crun.year, crun[k] * sc, c=COL['crunveg'], label='PICAL_crunveg (forest)')
        a.axvline(2100, c=MUTED, lw=1)
        if obs is not None:
            a.axhline(obs, c=INK, lw=1, ls=':')
            a.text(2041, obs, ' GIOMAS 1989-99', color=INK, fontsize=8, va='bottom')
        a.set(xlabel='model year', title=lab)
    ax[0].legend(fontsize=8, frameon=False)
    fig.suptitle('Production line across the 2100 switch, and the forested branch on the same restarts',
                 color=INK)
    fig.tight_layout()
    fig.savefig('plots/lpjg_canopy_climate_ice.png', dpi=140)


if __name__ == '__main__':
    veg()
    climate()
