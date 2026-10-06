"""Figures for report rounds 47-52 (2026-10-04 to 06).

  plots/kpp_gm1000_collapse.png   Antarctic September extent, net TOA and the September mixed-layer
                                  area south of 55S, 2170-2189/2209: TKE-only line, IDEMIX arm, KPP+GM1000
  plots/so_cap_mechanism.png      (a) salt and temperature parts of the September density step in the
                                  newly convecting region, arm against TKE-only; (b) September cap at the
                                  end of every run against PHC3
  plots/so_westerlies_u10.png     zonal-mean 10 m zonal wind over the ocean, model against ERA5

Usage: rounds47_52_figures.py     (reval environment, from the repo root)
"""
import os, glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
B = f'{P}/runtime/awiesm3-v3.4'
C = {'TKE only': '#1f4e79', 'TKE + IDEMIX': '#c0504d', 'KPP + GM 1000': '#e69f00', 'obs': 'k'}
RUNS = {'TKE only': ('tke_albsn082', 'PICAL_crunveg_tke_albsn082'), 'TKE + IDEMIX': ('idemix_albsn082', 'PICAL_crunveg_idemix_albsn082'),
        'KPP + GM 1000': ('kpp_gm1000', 'PICAL_crunveg_kpp_gm1000')}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})


def gmean(root, var, yr):
    with xr.open_dataset(f'{root}/outdata/oifs/atm_remapped_1m_{var}_{yr}-{yr}.nc') as ds:
        x = ds[var].mean([d for d in ds[var].dims if d not in ('lat', 'lon')]).mean('lon')
        w = np.cos(np.deg2rad(ds['lat']))
        return float((x * w).sum() / w.sum())


# ---------------------------------------------------------------- collapse
md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
lat = md['lat'].values; area = md['nod_area'].values[0]; S55 = lat < -55
fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
for lab, (tag, exp) in RUNS.items():
    s = pd.read_csv(f'{REPO}/data/clim/series_{tag}.csv').set_index('year').loc[2170:2209]
    ax[0].plot(s.index, s['siextents_m09'], color=C[lab], label=lab, lw=1.6)
    yrs = [y for y in range(2171, 2210) if os.path.exists(f'{B}/{exp}/outdata/oifs/atm_remapped_1m_tsr_{y}-{y}.nc')]
    toa = [(gmean(f'{B}/{exp}', 'tsr', y) + gmean(f'{B}/{exp}', 'ttr', y)) / 3600 for y in yrs]
    ax[1].plot(yrs, pd.Series(toa).rolling(3, center=True, min_periods=1).mean(), color=C[lab], lw=1.6)
    ym = [y for y in range(2171, 2210, 2) if os.path.exists(f'{B}/{exp}/outdata/fesom/MLD2.fesom.{y}.nc')]
    a = [area[S55 & (np.abs(xr.open_dataset(f'{B}/{exp}/outdata/fesom/MLD2.fesom.{y}.nc')['MLD2'].isel(time=8).values) > 1000)].sum() / 1e12 for y in ym]
    ax[2].plot(ym, a, color=C[lab], lw=1.6, marker='o', ms=3)
ax[0].set_ylabel('10$^6$ km$^2$'); ax[0].set_title('(a) Antarctic sea-ice extent, September', loc='left', fontsize=9); ax[0].legend(frameon=False, fontsize=8)
ax[1].axhline(0, color='0.6', lw=0.6); ax[1].set_ylabel('W m$^{-2}$'); ax[1].set_title('(b) net TOA (3-year running mean)', loc='left', fontsize=9)
ax[2].set_ylabel('10$^6$ km$^2$'); ax[2].set_title('(c) area with Sept. MLD > 1000 m, south of 55S', loc='left', fontsize=9)
for a_ in ax: a_.set_xlim(2170, 2209)
fig.tight_layout(); fig.savefig(f'{REPO}/plots/kpp_gm1000_collapse.png', dpi=200); plt.close(fig)

# ---------------------------------------------------------------- cap mechanism
D = pd.read_csv(f'{REPO}/data/clim/so_destabilisation_monthly.csv')
q = D[(D.region == 'NEW') & (D.month == 9)]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), gridspec_kw={'width_ratios': [1, 1.25]})
for run, lab, ls in (('tke', 'TKE only', '--'), ('arm', 'KPP + GM 1000', '-')):
    r = q[q.run == run].sort_values('year')
    ax[0].plot(r.year, r.dsig_S, color='#1f77b4', ls=ls, lw=1.6, label=f'salt part, {lab}')
    ax[0].plot(r.year, r.dsig_T, color='#d62728', ls=ls, lw=1.6, label=f'temperature part, {lab}')
    ax[0].plot(r.year, r.dsig, color='k', ls=ls, lw=1.6, label=f'net, {lab}')
ax[0].axhline(0, color='0.6', lw=0.6); ax[0].axhline(0.32, color='k', lw=0.8, ls=':'); ax[0].text(2176.6, 0.33, 'PHC3 net, 0.32', fontsize=7.5)
ax[0].set_ylabel('kg m$^{-3}$'); ax[0].legend(frameon=False, fontsize=7, ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.12))
ax[0].set_title('(a) September density step, 200-300 m minus 0-50 m', loc='left', fontsize=9)
R = pd.read_csv(f'{REPO}/data/clim/so_cap_drift_by_run.csv'); R = R[R.region == 'NEW']
order = ['gmramp 1200s', 'kpp_gm1000 KPP GM1000', 'gmferr1800 Ferreira', 'gmrd1800 Rossby radius', 'gm1500 1200s', 'crunveg 1200s GM2500 TKE+IDEMIX',
         'gm2500 1200s (ob1200)', 'idemix_albsn082 TKE+IDEMIX', 'ob1800 GM2500', 'tke_albsn082 TKE', 'gmhemi1800 TKE+IDEMIX', 'gmhemi_tke TKE']
names = ['GM off on fine cells', 'KPP + GM 1000', 'Ferreira', 'Rossby radius', 'GM 1500', 'GM 2500, 1200 s', 'GM 2500 (ob1200)', 'hemi, IDEMIX arm',
         'GM 2500, 1800 s', 'hemi, TKE only', 'hemi, TKE+IDEMIX', 'hemi, TKE 2141-59']
v = [R[R.run == o].sort_values('year').iloc[-3:]['step'].mean() for o in order]
col = ['#e69f00' if n in ('KPP + GM 1000',) else ('#7f7f7f' if i < 5 else '#1f4e79') for i, n in enumerate(names)]
ax[1].barh(range(len(v)), v, color=col); ax[1].set_yticks(range(len(v))); ax[1].set_yticklabels(names, fontsize=7.5)
ax[1].axvline(0.32, color='k', lw=1, ls=':'); ax[1].text(0.322, len(v) - 0.6, 'PHC3', fontsize=7.5)
ax[1].set_xlabel('kg m$^{-3}$'); ax[1].set_title('(b) September cap at the end of each run', loc='left', fontsize=9)
fig.tight_layout(); fig.savefig(f'{REPO}/plots/so_cap_mechanism.png', dpi=200); plt.close(fig)

# ---------------------------------------------------------------- westerlies
RUN = f'{B}/PICAL_crunveg_tke_albsn082/outdata/oifs'
fs = sorted(f for f in glob.glob(f'{RUN}/atm_remapped_1m_10u_*.nc') if 2195 <= int(f[-12:-8]) <= 2209)
m = xr.open_mfdataset(fs, combine='by_coords', use_cftime=True)['10u']
t = 'time_counter' if 'time_counter' in m.dims else 'time'
mc = m.mean(t).compute().sortby('lat')
l = xr.open_dataset(f'{RUN}/atm_remapped_1m_lsm_2209-2209.nc')['lsm']; l = l.isel({d: 0 for d in l.dims if d not in ('lat', 'lon')}).sortby('lat')
era = xr.open_dataset('/albedo/work/user/jstreffi/obs/era5/netcdf/u10.nc'); ev = [v_ for v_ in era.data_vars if era[v_].ndim >= 3][0]
e = era[ev].rename({k: v_ for k, v_ in (('latitude', 'lat'), ('longitude', 'lon')) if k in era[ev].dims})
ec = e.mean([d for d in e.dims if d not in ('lat', 'lon')]).sortby('lat').interp(lat=mc['lat'], lon=mc['lon']).compute()
oc = l < 0.5
fig, ax = plt.subplots(figsize=(5.2, 3.2))
ax.plot(mc['lat'], mc.where(oc).mean('lon'), color=C['TKE only'], lw=1.8, label='model, TKE-only line 2195-2209')
ax.plot(mc['lat'], ec.where(oc).mean('lon'), color='k', lw=1.4, label='ERA5 1989-2014')
ax.axhline(0, color='0.6', lw=0.6); ax.axvspan(-65, -50, color='0.92', zorder=0)
ax.set_xlim(-78, -20); ax.set_ylim(-5, 9); ax.set_xlabel('latitude'); ax.set_ylabel('m s$^{-1}$'); ax.legend(frameon=False, fontsize=8)
ax.set_title('Zonal-mean 10 m zonal wind over the ocean', loc='left', fontsize=9)
fig.tight_layout(); fig.savefig(f'{REPO}/plots/so_westerlies_u10.png', dpi=200); plt.close(fig)
print('written')
