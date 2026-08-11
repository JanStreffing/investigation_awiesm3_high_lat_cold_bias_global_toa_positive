"""TOA cloud forcing and surface shortwave bias, six arms, Southern Ocean to the tropics.

WHY THIS FRAME.  The 40-70S strip showed the levers acting as a near-uniform multiplier
rather than a spatial correction, and the global frame showed the tropics moving as much
as the Southern Ocean.  Neither view lets you see both regions at a readable scale at
once.  This one runs 80S to 30N: the whole scored SO band, the subtropical stratocumulus
decks that dominate the error field, and enough of the tropics to judge the guardrail.

WHY BOTH COLUMNS.  The campaign quotes the SO band mean from TOA SW CRE and the "SW RMSE"
from SURFACE net SW over ocean (eval_round10_A.py:190).  Those are different fields and
they disagree about which arm is best -- P5 has the surface nearly right and TOA 5 W/m2
short, LX3 the reverse.  Showing one without the other is how that went unnoticed.

READING THE STATS.  Each panel carries the SO band mean and the tropical mean, because a
lever that fixes the first by moving the second is a global dimming knob, not a Southern
Ocean fix.  The tropics sit BELOW CERES, so for them red (too much absorbed) is the error
and blue is worse, not better -- the opposite of the Southern Ocean reading.

DESIGN.  Diverging data about zero -> two hues, neutral midpoint, symmetric limits, one
shared scale down each column so panels are comparable by eye.  Per-panel autoscaling
would make the worst arm look the calmest.  Coastlines from the model's own land-sea mask;
cartopy is unusable on these nodes (its shapely wants a newer libstdc++).
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
OUT = ('/work/ab0246/a270092/postprocessing/'
       'investigation_awiesm3_high_lat_cold_bias_global_toa_positive/plots')
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/'
        'Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/outdata/oifs/'
        'atm_remapped_1m_lsm_1350-1350.nc')
ACC, Y0, Y1 = 3600.0, 1872, 1915
ARMS = [('(a) control', 'amip_pi_base'),
        ('(b) P5   adopted stack', 'amip_P5_swemin15'),
        ('(c) S4   + INPPMIN 50k', 'amip_S4_inppmin50000'),
        ('(d) LX1  + RSNOWLIN2', 'amip_LX1_long'),
        ('(e) LX3  + DMS + INPPMIN', 'amip_LX3_long'),
        ('(f) LY2  + overlap 0.10', 'amip_LY2_long')]
SO, TROP = (-65.0, -45.0), (-30.0, 30.0)
VIEW = (-80.0, 30.0)
LIM = 25.0

print(__doc__)
print('=' * 104)


def field(run, kind):
    need = ['tsr', 'tsrc'] if kind == 'cre' else ['ssr']
    acc, lat, lon = [], None, None
    for y in range(Y0, Y1 + 1):
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


with xr.open_dataset(CERESF) as c:
    OBS = {'cre': c['toa_cre_sw_clim'].values.mean(axis=0),
           'sfc': c['sfc_net_sw_all_clim'].values.mean(axis=0)}
    clat, clon = c['lat'].values, c['lon'].values
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values)
    if lsm.ndim == 3:
        lsm = lsm[0]
    llat, llon = d['lat'].values, d['lon'].values

COLS = [('TOA SW CRE', 'cre', False), ('surface net SW  (ocean only)', 'sfc', True)]
n = len(ARMS)
fig, axes = plt.subplots(n, 2, figsize=(17.0, 2.55 * n + 1.6), sharex=True, sharey=True)
cmap = plt.get_cmap('RdBu_r')
norm = matplotlib.colors.TwoSlopeNorm(vmin=-LIM, vcenter=0.0, vmax=LIM)

print(f'  {"arm":26s} {"field":14s} {"SO mean":>9s} {"SO RMSE":>9s} {"trop mean":>10s}')
for irow, (lab, run) in enumerate(ARMS):
    for jcol, (cname, kind, ocean) in enumerate(COLS):
        m, lat, lon = field(run, kind)
        ax = axes[irow, jcol]
        if m is None:
            ax.text(.5, .5, 'no data', ha='center', transform=ax.transAxes); continue
        ii = np.abs(clat[None, :] - lat[:, None]).argmin(axis=1)
        jj = np.abs(clon[None, :] - lon[:, None]).argmin(axis=1)
        li = np.abs(llat[None, :] - lat[:, None]).argmin(axis=1)
        lj = np.abs(llon[None, :] - lon[:, None]).argmin(axis=1)
        bias = m - OBS[kind][np.ix_(ii, jj)]
        msk = (lsm[np.ix_(li, lj)] <= 0.5) if ocean else np.ones_like(bias, bool)
        W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], bias.shape)
        s = ((lat >= SO[0]) & (lat < SO[1]))[:, None] & msk
        t = ((lat >= TROP[0]) & (lat < TROP[1]))[:, None] & msk
        som = float(np.average(bias[s], weights=W[s]))
        sor = float(np.sqrt(np.average(bias[s] ** 2, weights=W[s])))
        trm = float(np.average(bias[t], weights=W[t]))
        print(f'  {lab:26s} {cname[:12]:14s} {som:9.2f} {sor:9.2f} {trm:10.2f}')
        show = np.where(msk, bias, np.nan)
        ys = (lat >= VIEW[0]) & (lat <= VIEW[1])
        im = ax.pcolormesh(lon, lat[ys], show[ys], cmap=cmap, norm=norm,
                           shading='auto', rasterized=True)
        ax.contour(llon, llat, lsm, levels=[0.5], colors='k', linewidths=0.5)
        for la in (SO[0], SO[1]):
            ax.axhline(la, color='k', lw=0.9, ls='--', alpha=0.75)
        for la in TROP:
            ax.axhline(la, color='0.35', lw=0.7, ls=':', alpha=0.8)
        ax.set_ylim(*VIEW); ax.set_yticks([-75, -60, -45, -30, 0, 30])
        ax.tick_params(labelsize=8)
        ax.set_title(f'{lab}   |   SO {som:+.2f} (rmse {sor:.2f})   tropics {trm:+.2f}',
                     fontsize=9.5, pad=3, loc='left')
        if jcol == 0:
            ax.set_ylabel('lat', fontsize=8.5)
        if irow == 0:
            ax.text(0.5, 1.33, cname, transform=ax.transAxes, ha='center',
                    fontsize=12.5, fontweight='bold')
for ax in axes[-1]:
    ax.set_xlabel('longitude', fontsize=8.5)
    ax.set_xlim(0, 360); ax.set_xticks(np.arange(0, 361, 60))
fig.subplots_adjust(left=0.045, right=0.99, top=0.945, bottom=0.075,
                    hspace=0.34, wspace=0.05)
cax = fig.add_axes([0.34, 0.028, 0.32, 0.010])
cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cb.set_label('bias vs CERES  [W m$^{-2}$]     red = too much absorbed / too little '
             'reflection.   Dashed = 45-65S scored band;  dotted = 30S/30N.', fontsize=8.5)
cb.ax.tick_params(labelsize=8)
os.makedirs(OUT, exist_ok=True)
p = f'{OUT}/pattern_bias_multi.png'
fig.savefig(p, dpi=145, bbox_inches='tight')
print(f'\n  written: {p}')
