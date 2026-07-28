"""Double-check: AMIP vs CERES and CRUNCEP3 vs CERES, all on ONE grid with ONE land mask.

The claim under test is that CRUNCEP3 delivers ~22 W/m2 more surface net SW over
Siberian land in JJA than CERES EBAF does, which would make CRUNCEP3 -- not the model
-- the outlier, and would shrink the model's summer deficit from 33 to 11 W/m2.

That claim was originally assembled from two separate scripts using DIFFERENT grids and
DIFFERENT land masks (CRUNCEP on the TCO95 cell grid with a GRIB land-sea mask; CERES
interpolated to the model's regular grid with the coupled run's mask). A 22 W/m2
difference can easily be manufactured that way, so everything is redone here on the
CERES 1x1 grid with a single land mask derived from the model LSM.

Cross-checks built in:
  * CERES JJA from the monthly climatology file vs the dedicated JJA-mean file
  * sensitivity to the land-mask threshold (0.5 vs 0.8)
  * all three sources reported with identical area weighting
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

REPO = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/'
OBS = '/work/ab0246/a270092/obs/CERES/'
AMIP = '/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_base/outdata/oifs'
CRU = '/work/ab0995/a270270/input/cruncep_v7/CRUNCEP_noLPJG_1d_1901-1910_TCO95_calibrated_v3.nc'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
ACC = 3600.0
JJA = [5, 6, 7]
BOXES = {'Siberia': (55, 75, 60, 180), 'E. Siberia': (55, 75, 90, 160)}

# ---- target grid: CERES 1x1 -------------------------------------------------
ce = xr.open_dataset(OBS + 'CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc')
LAT, LON = ce.lat.values, ce.lon.values                     # 0.5..359.5
ceres = ce['sfc_net_sw_all_clim'].values                    # (12, 180, 360)
ceres_clr = ce['sfc_net_sw_clr_t_clim'].values
ce.close()

# ---- AMIP: regular lat-lon -> CERES grid ------------------------------------
acc = []
for y in range(1872, 1880):
    ds = xr.open_dataset(f'{AMIP}/atm_remapped_1m_ssr_1m_{y}-{y}.nc')
    acc.append(ds['ssr'].values / ACC)
    mlat, mlon = ds['ssr'].lat.values, ds['ssr'].lon.values
    ds.close()
amip_native = np.mean(acc, axis=0)
amip = xr.DataArray(amip_native, dims=('t', 'lat', 'lon'),
                    coords={'t': np.arange(12), 'lat': mlat, 'lon': mlon}
                    ).interp(lat=LAT, lon=LON, kwargs={'fill_value': None}).values

# ---- model LSM -> CERES grid (the single land mask used for everything) ------
lsm_native = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0)
LSM = lsm_native.interp(lat=LAT, lon=LON, kwargs={'fill_value': None}).values

# ---- CRUNCEP3: unstructured cells -> CERES grid, spherical nearest neighbour --
cr = xr.open_dataset(CRU)
clat, clon = cr.lat.values, cr.lon.values
mon = cr['time_counter'].dt.month.values
cru_mon = np.stack([cr['rsns'].values[mon == m].mean(axis=0) for m in range(1, 13)])
cr.close()


def xyz(la, lo):
    la, lo = np.deg2rad(la), np.deg2rad(lo)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


tree = cKDTree(xyz(clat, clon))
G_LON, G_LAT = np.meshgrid(LON, LAT)
_, idx = tree.query(xyz(G_LAT.ravel(), G_LON.ravel()))
cru = cru_mon[:, idx].reshape(12, len(LAT), len(LON))

# ---- area-weighted box means ------------------------------------------------
def boxmean(f, box, months, thr=0.5):
    la0, la1, lo0, lo1 = box
    lon180 = ((LON + 180) % 360) - 180
    yi, xi = (LAT >= la0) & (LAT <= la1), (lon180 >= lo0) & (lon180 <= lo1)
    sub = f[np.ix_(months, np.where(yi)[0], np.where(xi)[0])]
    m = LSM[np.ix_(np.where(yi)[0], np.where(xi)[0])] > thr
    w = np.where(m, np.cos(np.deg2rad(LAT[yi]))[:, None] * np.ones(m.shape), 0.0)
    return float((sub * w).sum() / (w.sum() * len(months)))


print("Surface NET shortwave, JJA, land only, common CERES 1x1 grid + single LSM\n")
print(f"  {'box':12s} {'CERES':>8s} {'AMIP':>8s} {'CRUNCEP3':>9s} | {'AMIP-CER':>9s} {'CRU-CER':>9s}")
for name, box in BOXES.items():
    c, a, u = (boxmean(x, box, JJA) for x in (ceres, amip, cru))
    print(f"  {name:12s} {c:8.1f} {a:8.1f} {u:9.1f} | {a-c:+9.1f} {u-c:+9.1f}")

print("\n  land-mask sensitivity (Siberia, threshold 0.5 -> 0.8):")
for thr in (0.5, 0.8):
    c, a, u = (boxmean(x, BOXES['Siberia'], JJA, thr) for x in (ceres, amip, cru))
    print(f"    thr={thr}: CERES {c:6.1f}  AMIP {a:6.1f} ({a-c:+5.1f})  CRUNCEP {u:6.1f} ({u-c:+5.1f})")

# cross-check the CERES JJA value against the dedicated JJA-mean file
jf = OBS + 'CERES_EBAF_Ed4.1_Subset_200003-202106_JJA.nc'
if os.path.exists(jf):
    dj = xr.open_dataset(jf)
    v = [k for k in dj.data_vars if 'sfc_net_sw_all' in k]
    if v:
        arr = dj[v[0]].values
        arr = arr[None] if arr.ndim == 2 else arr.mean(axis=0)[None]
        print(f"\n  cross-check, CERES JJA from dedicated JJA file: "
              f"{boxmean(arr, BOXES['Siberia'], [0]):.1f} "
              f"(vs {boxmean(ceres, BOXES['Siberia'], JJA):.1f} from the monthly climatology)")
    dj.close()

# ---- figure -----------------------------------------------------------------
SURF, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.linewidth': 0.6,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.labelcolor': INK2, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF})
yi = (LAT >= 40) & (LAT <= 80)
xi = (LON >= 0) & (LON <= 190)
ext = [LON[xi][0], LON[xi][-1], LAT[yi][0], LAT[yi][-1]]
land = LSM[np.ix_(np.where(yi)[0], np.where(xi)[0])] > 0.5
cj = ceres[JJA][:, yi][:, :, xi].mean(0)
aj = amip[JJA][:, yi][:, :, xi].mean(0)
uj = cru[JJA][:, yi][:, :, xi].mean(0)

fig = plt.figure(figsize=(12.4, 6.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.95], hspace=0.38, wspace=0.22)
for k, (f, ttl, kw) in enumerate([
        (np.where(land, cj, np.nan), 'CERES EBAF, JJA surface net SW', dict(cmap='YlOrBr', vmin=120, vmax=210)),
        (np.where(land, aj - cj, np.nan), 'AMIP $-$ CERES', dict(cmap='RdBu_r', vmin=-40, vmax=40)),
        (np.where(land, uj - cj, np.nan), 'CRUNCEP3 $-$ CERES', dict(cmap='RdBu_r', vmin=-40, vmax=40))]):
    ax = fig.add_subplot(gs[0, k])
    im = ax.imshow(f, origin='lower', extent=ext, aspect='auto', **kw)
    ax.set_title(ttl, fontsize=9.5, color=INK, loc='left', fontweight='bold', pad=5)
    ax.set_xlabel('lon'); ax.set_ylabel('lat' if k == 0 else '')
    for (la0, la1, lo0, lo1), c in [((55, 75, 60, 180), '#0b0b0b')]:
        ax.plot([lo0, lo1, lo1, lo0, lo0], [la0, la0, la1, la1, la0], color=c, lw=1.2, ls='--')
    fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02).ax.tick_params(labelsize=7)

ax = fig.add_subplot(gs[1, :2])
mm = np.arange(1, 13)
for f, lab, c in [(ceres, 'CERES EBAF', '#2a78d6'), (amip, 'AMIP', '#eb6834'),
                  (cru, 'CRUNCEP3', '#1baf7a')]:
    ax.plot(mm, [boxmean(f, BOXES['Siberia'], [i]) for i in range(12)],
            color=c, lw=2.2, marker='o', ms=3.5, mec=SURF, mew=0.8, label=lab)
ax.axvspan(5.5, 8.5, color=MUTED, alpha=0.10, lw=0)
ax.set_xticks(mm); ax.set_xticklabels(list('JFMAMJJASOND'))
ax.set_ylabel('surface net SW [W m$^{-2}$]')
ax.set_title('Siberia land (55--75N, 60--180E), seasonal cycle', fontsize=9.5, color=INK,
             loc='left', fontweight='bold', pad=5)
ax.grid(axis='y', color=MUTED, alpha=0.22, lw=0.5); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=8.4, labelcolor=INK2)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)

ax = fig.add_subplot(gs[1, 2])
for f, lab, c in [(amip, 'AMIP $-$ CERES', '#eb6834'), (cru, 'CRUNCEP3 $-$ CERES', '#1baf7a')]:
    ax.plot(mm, [boxmean(f, BOXES['Siberia'], [i]) - boxmean(ceres, BOXES['Siberia'], [i])
                 for i in range(12)], color=c, lw=2.2, marker='o', ms=3.5, mec=SURF, mew=0.8, label=lab)
ax.axhline(0, color=MUTED, lw=0.8); ax.axvspan(5.5, 8.5, color=MUTED, alpha=0.10, lw=0)
ax.set_xticks(mm); ax.set_xticklabels(list('JFMAMJJASOND'))
ax.set_ylabel('difference [W m$^{-2}$]')
ax.set_title('difference from CERES', fontsize=9.5, color=INK, loc='left', fontweight='bold', pad=5)
ax.grid(axis='y', color=MUTED, alpha=0.22, lw=0.5); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=8.0, labelcolor=INK2)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)

fig.suptitle('Verification: AMIP and CRUNCEP3 against CERES EBAF on a common grid',
             x=0.008, ha='left', fontsize=11.5, fontweight='bold', color=INK, y=0.985)
fig.text(0.008, 0.938, 'Surface net shortwave, land only, all three regridded to CERES 1x1 with a single '
         'land mask. Dashed box = Siberia (55--75N, 60--180E).', ha='left', fontsize=8.2, color=INK2)
fig.savefig(REPO + 'plots/verify_ceres_amip_cruncep.png', dpi=170, bbox_inches='tight')
print('\nwrote plots/verify_ceres_amip_cruncep.png')
