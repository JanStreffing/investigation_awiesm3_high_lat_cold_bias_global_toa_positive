"""GLOBAL view of the two candidate configurations, TOA cloud forcing and surface shortwave.

WHY GLOBAL.  The Southern Ocean strip showed the levers acting as a near-uniform multiplier
rather than a spatial correction, and showed the TOA and SURFACE constraints disagreeing --
P5 has the surface about right and TOA 5 W/m2 short, LX3 the reverse.  Neither of those is
a Southern Ocean fact if the levers are also acting elsewhere, and RSNOWLIN2 in particular
is a TEMPERATURE-dependent ice-cloud knob, so it must act wherever cold cloud exists.  This
frame shows where else the energy moved.

THE PARADOX TO RESOLVE.  Against CERES, the 45-65S band-mean SW CRE error goes from
+7.85 W/m2 (control) to essentially closed in LX3 -- and the spatial RMSE goes the WRONG
way: control 6.877, P5 4.938, LX3 7.329.  LX3 is worse than doing nothing on pattern while
being much better on the mean.  That can only happen if the extra reflection lands in the
wrong places, so the question is a map question.

WHAT IS PLOTTED, and the resolution turns out to be that TWO DIFFERENT FIELDS were being
compared.  The band mean quoted in the campaign is TOA SW CRE; the "SW RMSE" that degrades
is SURFACE net shortwave over ocean (eval_round10_A.py:190, rmse(ssrm, csw, ...) with
ssrm = model ssr and csw = CERES sfc_net_sw_all).  So the left column is TOA SW CRE bias
and the right column is surface net SW bias, same arms, same scale within each column.
LX3 improves the TOA CRE pattern monotonically (RMSE 11.45 -> 8.78 -> 7.74) while the
surface field is what worsens.  Plotting them together is the only way to see that.

DESIGN, and why each choice is forced:
  * DIVERGING data (an error about zero) -> two hues with a NEUTRAL midpoint, symmetric
    limits.  Never a rainbow: rainbow maps invent structure at their hue boundaries and
    hide it inside bands, which is exactly the thing being judged here.
  * ONE shared colour scale and ONE symmetric limit across all three maps.  Per-panel
    autoscaling would make the worst arm look the calmest.
  * The band mean and the spatial RMSE are printed ON each panel, because the whole point
    is that those two numbers disagree.
  * Panel (d) is the argument: if LX3's histogram is centred better but broader than P5's,
    the mean improved by cancellation of larger errors of both signs.

PERIOD.  1872-1915, the campaign standard, against the CERES EBAF climatology.  Both are
climatologies of different epochs; the epoch offset on this term was measured at
-0.07 W/m2 (round 22) and is negligible next to the effects here.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# NO CARTOPY: its shapely needs a GLIBCXX newer than this node's libstdc++, so it
# imports and dies.  A plain longitude-latitude strip is arguably better here anyway --
# the Southern Ocean band is zonal, and a lon-lat plot shows the longitude structure
# (the sectors) directly instead of wrapping it around a pole.  Coastlines are contoured
# from the model's OWN land-sea mask, which needs no external geometry library.

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
OUT = ('/work/ab0246/a270092/postprocessing/'
       'investigation_awiesm3_high_lat_cold_bias_global_toa_positive/report/plots')
ACC, Y0, Y1 = 3600.0, 1872, 1915
ARMS = [('P5  adopted stack', 'amip_P5_swemin15'),
        ('LX3  +rsnow+DMS+INP', 'amip_LX3_long')]
BAND = (-65.0, -45.0)          # the scored band
VIEW = (-90.0, 90.0)           # GLOBAL: the SO levers are not confined to the SO
LIM = 30.0                     # symmetric colour limit, W/m2
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/'
        'Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/outdata/oifs/'
        'atm_remapped_1m_lsm_1350-1350.nc')

print(__doc__)
print('=' * 96)


def model_field(run, kind):
    """kind='cre' -> TOA SW CRE (tsr-tsrc); kind='sfc' -> surface net SW (ssr)."""
    acc, lat, lon = [], None, None
    for y in range(Y0, Y1 + 1):
        if kind == 'cre':
            need = ['tsr', 'tsrc']
        else:
            need = ['ssr']
        fs = [f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc' for v in need]
        if not all(os.path.exists(f) for f in fs):
            continue
        vals = []
        for v, f in zip(need, fs):
            with xr.open_dataset(f, decode_times=False) as d:
                vals.append(d[v].values / ACC)
                lat, lon = d['lat'].values, d['lon'].values
        a = vals[0] - vals[1] if kind == 'cre' else vals[0]
        if a.shape[0] == 12:
            acc.append(a)
    return (np.mean(acc, axis=0).mean(axis=0) if acc else None), lat, lon


def model_swcre(run):
    acc, lat, lon = [], None, None
    for y in range(Y0, Y1 + 1):
        f1 = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_tsr_1m_{y}-{y}.nc'
        f2 = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_tsrc_1m_{y}-{y}.nc'
        if not (os.path.exists(f1) and os.path.exists(f2)):
            continue
        with xr.open_dataset(f1, decode_times=False) as d:
            a = d['tsr'].values / ACC
            lat, lon = d['lat'].values, d['lon'].values
        with xr.open_dataset(f2, decode_times=False) as d:
            b = d['tsrc'].values / ACC
        if a.shape[0] == 12:
            acc.append(a - b)
    return (np.mean(acc, axis=0).mean(axis=0) if acc else None), lat, lon


with xr.open_dataset(CERESF) as c:
    ce_cre = c['toa_cre_sw_clim'].values.mean(axis=0)
    ce_sfc = c['sfc_net_sw_all_clim'].values.mean(axis=0)
    clat, clon = c['lat'].values, c['lon'].values

with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values)
    if lsm.ndim == 3:
        lsm = lsm[0]
    llat, llon = d['lat'].values, d['lon'].values

COLS = [('TOA SW CRE', 'cre', ce_cre), ('surface net SW', 'sfc', ce_sfc)]
LIMG = 20.0
fig, axes = plt.subplots(2, 2, figsize=(17.5, 8.6), sharex=True, sharey=True)
cmap = plt.get_cmap('RdBu_r')
norm = matplotlib.colors.TwoSlopeNorm(vmin=-LIMG, vcenter=0.0, vmax=LIMG)

print(f'  {"arm":22s} {"field":14s} {"global mean":>12s} {"global RMSE":>12s} '
      f'{"45-65S mean":>12s} {"tropics mean":>13s}')
for irow, (lab, run) in enumerate(ARMS):
    for jcol, (cname, kind, obs_src) in enumerate(COLS):
        m, lat, lon = model_field(run, kind)
        ax = axes[irow, jcol]
        if m is None:
            ax.text(.5, .5, 'no data', ha='center', transform=ax.transAxes); continue
        ii = np.abs(clat[None, :] - lat[:, None]).argmin(axis=1)
        jj = np.abs(clon[None, :] - lon[:, None]).argmin(axis=1)
        bias = m - obs_src[np.ix_(ii, jj)]
        W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], bias.shape)
        gm = float(np.average(bias, weights=W))
        gr = float(np.sqrt(np.average(bias**2, weights=W)))
        so = (lat >= -65) & (lat < -45); tr = (lat >= -30) & (lat < 30)
        som = float(np.average(bias[so], weights=W[so]))
        trm = float(np.average(bias[tr], weights=W[tr]))
        print(f'  {lab:22s} {cname:14s} {gm:12.3f} {gr:12.3f} {som:12.3f} {trm:13.3f}')
        im = ax.pcolormesh(lon, lat, bias, cmap=cmap, norm=norm, shading='auto',
                           rasterized=True)
        ax.contour(llon, llat, lsm, levels=[0.5], colors='k', linewidths=0.5)
        for la in (-65, -45, -30, 30):
            ax.axhline(la, color='0.25', lw=0.7, ls='--', alpha=0.7)
        ax.set_ylim(-90, 90); ax.set_yticks(np.arange(-90, 91, 30))
        ax.tick_params(labelsize=8.5)
        ax.set_title(f'{lab}  |  {cname}    global mean {gm:+.2f}   RMSE {gr:.2f}   '
                     f'[45-65S {som:+.2f},  tropics {trm:+.2f}]', fontsize=10, pad=4,
                     loc='left')
        if jcol == 0:
            ax.set_ylabel('latitude', fontsize=9)
for ax in axes[-1]:
    ax.set_xlabel('longitude', fontsize=9); ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 60))
fig.subplots_adjust(left=0.045, right=0.985, top=0.955, bottom=0.13, hspace=0.20, wspace=0.06)
cax = fig.add_axes([0.33, 0.055, 0.36, 0.016])
cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cb.set_label('bias vs CERES  [W m$^{-2}$]   (red = too much absorbed / too little reflection)',
             fontsize=9)
cb.ax.tick_params(labelsize=8.5)
os.makedirs(OUT, exist_ok=True)
pth = f'{OUT}/global_pattern_bias.png'
fig.savefig(pth, dpi=150, bbox_inches='tight')
print(f'\n  written: {pth}')
