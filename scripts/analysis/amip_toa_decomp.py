"""Decompose the AMIP atmosphere-only TOA imbalance (+0.67 W/m2).

Where does the residual atmospheric radiative error live -- which hemisphere,
which latitude band, SW or LW, clear-sky or cloud? This determines which knob
can close it without undoing the boreal-summer fix (which needs LESS cloud,
i.e. pushes TOA the wrong way).

Sign conventions (IFS, all fluxes positive DOWNWARD at TOA):
  tsr  = TOA net solar   = incident - reflected     (positive)
  ttr  = TOA net thermal = -OLR                     (negative)
  tisr = TOA incident solar                         (positive)
  tsrc/ttrc = the same, clear-sky
  SW CRE = tsr - tsrc   (negative: clouds reflect)
  LW CRE = ttr - ttrc   (positive: clouds trap OLR)
  planetary albedo = 1 - tsr/tisr

Benchmark: CERES-EBAF global climatology
  absorbed SW 240.5, OLR 240.2, albedo 0.293,
  SW CRE -45.4, LW CRE +25.8, net CRE -19.6 W/m2
PI-vs-present differences in CRE are a few tenths, so CERES is a fair yardstick
at the several-W/m2 level for structure and magnitude.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ACC = 3600.0
REPO = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/'
AMIP = '/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_base/outdata/oifs'
C09C = '/work/bb1469/a270092/runtime/awiesm3-v3.4/Tuning_test_09C_06V_CRUNCEPinit_newSeaIce/outdata/oifs'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
VARS = ['tsr', 'ttr', 'tsrc', 'ttrc', 'tisr']

CERES = dict(aswabs=240.5, olr=240.2, alb=0.293, swcre=-45.4, lwcre=25.8, netcre=-19.6)
# Tiling bands: their weighted contributions must sum to the global mean.
BANDS = [('90S-60S', -90, -60), ('60S-30S', -60, -30), ('tropics 30S-30N', -30, 30),
         ('30N-60N', 30, 60), ('60N-90N', 60, 90)]
# Diagnostic overlay only -- OVERLAPS the tiling bands, so it is excluded from the sum.
OVERLAY = [('SO 65S-45S (overlay)', -65, -45)]


def load(dirpath, infix, y0, y1):
    """Annual-mean 2-D field per variable, averaged over the year range."""
    out, lat, lon = {}, None, None
    for v in VARS:
        acc = []
        for y in range(y0, y1 + 1):
            f = f'{dirpath}/atm_remapped_1m_{v}{infix}_{y}-{y}.nc'
            if not os.path.exists(f):
                continue
            ds = xr.open_dataset(f)
            da = ds[v].mean('time_counter') / ACC
            if lat is None:
                lat, lon = da.lat.values, da.lon.values
            acc.append(da.values)
            ds.close()
        out[v] = np.mean(acc, axis=0)
    return out, lat, lon


def wmean(field, lat, mask=None):
    w = np.cos(np.deg2rad(lat))[:, None] * np.ones(field.shape)
    if mask is not None:
        w = np.where(mask, w, 0.0)
    return float((field * w).sum() / w.sum())


def report(tag, d, lat, lsm):
    net = d['tsr'] + d['ttr']
    netc = d['tsrc'] + d['ttrc']
    swcre = d['tsr'] - d['tsrc']
    lwcre = d['ttr'] - d['ttrc']
    g = lambda f, m=None: wmean(f, lat, m)
    alb = 1 - g(d['tsr']) / g(d['tisr'])
    print(f"\n{'='*76}\n{tag}\n{'='*76}")
    print(f"  net TOA          {g(net):8.3f}   clear-sky net {g(netc):8.3f}   "
          f"=> cloud contributes {g(net)-g(netc):+7.3f}")
    print(f"  absorbed SW      {g(d['tsr']):8.3f}   (CERES {CERES['aswabs']:.1f}, "
          f"diff {g(d['tsr'])-CERES['aswabs']:+.2f})")
    print(f"  OLR              {-g(d['ttr']):8.3f}   (CERES {CERES['olr']:.1f}, "
          f"diff {-g(d['ttr'])-CERES['olr']:+.2f})")
    print(f"  planetary albedo {alb:8.4f}   (CERES {CERES['alb']:.3f}, diff {alb-CERES['alb']:+.4f})")
    print(f"  SW CRE           {g(swcre):8.3f}   (CERES {CERES['swcre']:.1f}, "
          f"diff {g(swcre)-CERES['swcre']:+.2f})")
    print(f"  LW CRE           {g(lwcre):8.3f}   (CERES {CERES['lwcre']:.1f}, "
          f"diff {g(lwcre)-CERES['lwcre']:+.2f})")
    print(f"  net CRE          {g(swcre)+g(lwcre):8.3f}   (CERES {CERES['netcre']:.1f}, "
          f"diff {g(swcre)+g(lwcre)-CERES['netcre']:+.2f})")
    land = lsm > 0.5
    print(f"  land / ocean net TOA:  {g(net, land):7.3f} / {g(net, ~land):7.3f}")

    # contribution of each band to the GLOBAL mean (area-weighted share)
    print(f"\n  {'band':22s} {'net TOA':>9s} {'contrib':>9s} {'SW CRE':>9s} {'LW CRE':>9s}")
    total = 0.0
    for name, a, b in BANDS:
        sel = (lat >= a) & (lat < b)
        m = np.zeros(net.shape, bool); m[sel, :] = True
        contrib = wmean(np.where(m, net, 0.0), lat)   # area-weighted share of the global mean
        total += contrib
        print(f"  {name:22s} {wmean(net,lat,m):9.3f} {contrib:9.3f} "
              f"{wmean(swcre,lat,m):9.2f} {wmean(lwcre,lat,m):9.2f}")
    print(f"  {'sum of tiling bands':22s} {'':9s} {total:9.3f}  <- must equal global net TOA")
    for name, a, b in OVERLAY:
        sel = (lat >= a) & (lat < b)
        m = np.zeros(net.shape, bool); m[sel, :] = True
        print(f"  {name:22s} {wmean(net,lat,m):9.3f} {'--':>9s} "
              f"{wmean(swcre,lat,m):9.2f} {wmean(lwcre,lat,m):9.2f}")

    # clear-sky components, for separating a cloud error from a clear-sky one
    olrc = -g(d["ttr"]) + g(lwcre)   # OLR_clear = OLR + LW CRE
    print(f"\n  clear-sky absorbed SW {g(d['tsr'])-g(swcre):8.2f}  (CERES ~287.5)")
    print(f"  clear-sky OLR         {olrc:8.2f}  (CERES ~266.5)")
    return net, swcre, lwcre


amip, lat, lon = load(AMIP, '_1m', 1872, 1879)
c09c, _, _ = load(C09C, '', 1370, 1379)
lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values

na, sa, la = report('AMIP 1872-79  (prescribed SST, prescribed vegetation)', amip, lat, lsm)
nc, sc, lc = report('09C 1370-79   (coupled, 06V + newSeaIce)', c09c, lat, lsm)

# ---------------- zonal figure ----------------
SURF, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.linewidth': 0.6,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.labelcolor': INK2, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF})
BLUE, ORANGE = '#2a78d6', '#eb6834'
fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.9))
for ax, (fa, fc, ttl) in zip(axes, [(na, nc, 'net TOA'), (sa, sc, 'SW cloud radiative effect'),
                                    (la, lc, 'LW cloud radiative effect')]):
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.plot(lat, fa.mean(axis=1), color=BLUE, lw=2.0, label='AMIP (atmosphere only)')
    ax.plot(lat, fc.mean(axis=1), color=ORANGE, lw=2.0, label='09C (coupled)')
    ax.set_title(ttl, fontsize=9.5, color=INK, loc='left', fontweight='bold', pad=6)
    ax.set_xlabel('latitude'); ax.set_xlim(-90, 90); ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.grid(axis='y', color=MUTED, alpha=0.22, lw=0.5); ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
axes[0].set_ylabel('W m$^{-2}$')
axes[0].legend(frameon=False, fontsize=8.2, labelcolor=INK2, loc='lower center')
fig.suptitle('Where the atmosphere-only TOA imbalance lives — zonal means',
             x=0.008, ha='left', fontsize=11.5, fontweight='bold', color=INK, y=1.02)
fig.text(0.008, 0.93, 'AMIP 1872-79 (prescribed SST) vs coupled 09C 1370-79. Positive net TOA = energy into '
         'the system.', ha='left', fontsize=8.2, color=INK2)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(REPO + 'plots/amip_toa_decomposition_zonal.png', dpi=170, bbox_inches='tight')
print('\nwrote plots/amip_toa_decomposition_zonal.png')
