"""Ocean mixing baseline: mixed layer depth, model vs WOA18, with the definitions forced to agree.

WHY.  The campaign's standing hypothesis is that the model removes heat from the surface
faster than it should: it absorbs a positive TOA imbalance while sitting cold at the
surface and filling the interior.  If that is right the mixed layer is too DEEP, and the
place to look before touching a single mixing parameter is the mixed layer depth itself.

THE DEFINITION PROBLEM, WHICH IS THE WHOLE POINT OF THIS SCRIPT.  "Mixed layer depth" is
not one quantity.  Four definitions are in play here and they do not agree:

  WOA18 M_an   potential density at 10 m as the reference; MLD is where sigma_theta
               exceeds that reference by 0.125 kg/m3.  (Stated verbatim in the file's
               `summary` attribute -- not inferred.)
  FESOM MLD2   same 0.125 kg/m3 threshold, but referenced to the SURFACE level
               (rhopot(nzmin), mid-depth 2.5 m on this mesh), not to 10 m.
               oce_ale_pressure_bv.F90:463.
  FESOM MLD1   a completely different criterion -- Large et al. (1997), the shallowest
               depth where the local buoyancy gradient matches the maximum gradient
               above it.  Not comparable to WOA at all.
  FESOM MLD3   MLD2's criterion with the CMOR 0.03 kg/m3 threshold.

Comparing FESOM MLD2 against WOA M_an therefore compares a surface-referenced MLD with a
10 m-referenced one.  Where a fresh or warm skin sits above 10 m the surface-referenced
version trips its threshold earlier and reads SHALLOWER -- i.e. it flatters the model on
exactly the hypothesis under test.  So this script does not use FESOM's diagnostic for the
comparison.  It recomputes the model MLD offline from the model's own monthly T and S with
WOA's criterion:

    sigma_theta(z) - sigma_theta(10 m) >= 0.125 kg/m3, linearly interpolated in z,

with sigma_theta from EOS-80 at 0 dbar reference (`seawater.dens0`), which is the equation
of state WOA18's sigma-theta is built on.  FESOM's MLD2 is carried alongside as a separate
series so the size of the definitional effect is a measured number rather than a worry.

WHAT STILL DOES NOT MATCH, and cannot be made to.
  * WOA computes MLD per PROFILE and then objectively analyses the result.  This computes
    MLD from a monthly-MEAN model profile.  MLD is a nonlinear functional of the profile,
    so in regions of high sub-monthly variability the model number is biased shallow
    relative to obs.  Mitigated (not removed) by computing MLD for each month of each year
    separately and averaging the MLDs, never the profiles.
  * WOA18 is 1981-2010.  The model arm is an 1850 control.
  * WOA has real coverage problems in the Southern Ocean in winter; M_dd (observation
    count) is carried through so the maps can be masked on it.

METHOD.  MLD is computed on the native unstructured mesh, in the mesh's own geometry, and
only then binned to WOA's 1 degree grid -- MLD is a nonlinear functional of the profile, so
regridding T and S first and computing afterwards would be a different quantity.  Zonal
means are additionally computed on each grid in its own geometry, with no regridding at
all, as a cross-check that the binning does not carry the result.

MESH GUARD.  11X and later run the corrected CORE3 (220509 nodes).  core3_beta (211567) is
a different mesh and mixing them silently is a known way to get a wrong answer here.

Usage:  ARM=11X YEARS=1380-1389 python3 scripts/analysis/mld_baseline_vs_woa18.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import glob
import sys
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')
import seawater
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ----------------------------------------------------------------------------- config
R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
WOA_DIR = '/work/ab0246/a270092/obs/WOA'
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PLOTS = os.path.join(REPO, 'plots')

ARM = os.environ.get('ARM', '11X')
_y = os.environ.get('YEARS', '1380-1389').split('-')
YEARS = list(range(int(_y[0]), int(_y[-1]) + 1))

# The corrected CORE3.  core3_beta is a DIFFERENT mesh (211567 nodes); see
# core3-mesh-swap-beta-vs-new.
MESH = os.environ.get('MESH', '/work/ab0246/a270092/input/fesom2/core3')
N_NODES_CORE3 = 220509

SIGMA_CRIT = 0.125      # kg/m3, WOA18's threshold
REF_DEPTH = 10.0        # m, WOA18's reference depth
# M_an is the OBJECTIVELY ANALYSED field and is complete over the ocean (41088 of 64800
# 1-degree cells).  M_dd is the raw profile count and is sparse -- in September, only 263
# of 4007 ocean cells south of 60S hold a single profile.  Masking the analysed field on
# M_dd therefore deletes most of the ocean, so it is NOT used as a mask; M_dd is reported
# per band instead, because "the analysis says 210 m" means something different where it
# rests on 6.6 % coverage than where it rests on 24 %.
MIN_OBS = 0

BANDS = [('90-60S  Antarctic/SO', -90, -60), ('60-45S  subantarctic', -60, -45),
         ('45-30S', -45, -30), ('30S-30N  tropics', -30, 30), ('30-45N', 30, 45),
         ('45-60N  subpolar NA', 45, 60), ('60-90N  Arctic/Nordic', 60, 90)]
BANDS = BANDS[::-1]          # report north to south
MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']


# ------------------------------------------------------------------------------- mesh
def load_mesh():
    with xr.open_dataset(f'{MESH}/mesh.nc', decode_times=False) as m:
        lon = m['lon'].values.astype('f8')
        lat = m['lat'].values.astype('f8')
        area = m['cell_area'].values.astype('f8')
        depth = m['depth'].values.astype('f8')          # mid-depths, positive down
    n = lon.size
    if n != N_NODES_CORE3:
        sys.exit(f'MESH GUARD: {MESH} has {n} nodes, expected {N_NODES_CORE3} (corrected '
                 f'CORE3). core3_beta is 211567 and is NOT interchangeable.')
    return lon, lat, area, depth


# --------------------------------------------------------------- MLD, WOA's definition
def mld_woa_criterion(temp, salt, depth):
    """MLD [m, positive down] from potential temperature and practical salinity.

    temp, salt : (nod2, nz) with NaN below the bottom
    depth      : (nz,) mid-depths, positive down, increasing
    """
    sig = seawater.dens0(salt, temp) - 1000.0            # EOS-80 sigma_theta at 0 dbar
    nz = depth.size

    # reference sigma at exactly REF_DEPTH, linearly interpolated between the bracketing
    # mid-depth levels (on this mesh: 7.5 m and 15 m).
    k = int(np.searchsorted(depth, REF_DEPTH))
    if depth[0] >= REF_DEPTH:
        ref = sig[:, 0].copy()
    else:
        k = min(max(k, 1), nz - 1)
        w = (REF_DEPTH - depth[k - 1]) / (depth[k] - depth[k - 1])
        ref = sig[:, k - 1] + w * (sig[:, k] - sig[:, k - 1])
        shallow = np.isnan(ref)                          # column shallower than 10 m
        if shallow.any():
            ref[shallow] = sig[shallow, 0]

    excess = sig - ref[:, None]                          # (nod2, nz)
    valid = ~np.isnan(sig)
    over = (excess >= SIGMA_CRIT) & valid

    nnod = sig.shape[0]
    mld = np.full(nnod, np.nan)
    hit = over.any(axis=1)
    first = np.argmax(over, axis=1)                      # first True; 0 where none

    idx = np.where(hit)[0]
    j = first[idx]
    j_lo = np.maximum(j - 1, 0)
    e_hi = excess[idx, j]
    e_lo = excess[idx, j_lo]
    z_hi = depth[j]
    z_lo = depth[j_lo]
    denom = e_hi - e_lo
    frac = np.where(np.abs(denom) > 1e-12, (SIGMA_CRIT - e_lo) / denom, 0.0)
    frac = np.clip(frac, 0.0, 1.0)
    mld[idx] = z_lo + frac * (z_hi - z_lo)
    mld[idx[j == 0]] = depth[0]                          # already over at the top level

    # never reached: mixed to the bottom
    nlev = valid.sum(axis=1)
    wet = nlev > 0
    bottom = np.full(nnod, np.nan)
    bottom[wet] = depth[np.clip(nlev[wet] - 1, 0, nz - 1)]
    to_bottom = wet & ~hit
    mld[to_bottom] = bottom[to_bottom]
    mld[~wet] = np.nan
    return mld, to_bottom, bottom


# ------------------------------------------------------------------------ model fields
def model_monthly_mld(root, depth):
    """12-month climatology of the model MLD under WOA's criterion, plus FESOM's MLD2."""
    acc = acc2 = accb = None
    nyr = 0
    for y in YEARS:
        ft = f'{root}/outdata/fesom/temp.fesom.{y}.nc'
        fs = f'{root}/outdata/fesom/salt.fesom.{y}.nc'
        if not (os.path.exists(ft) and os.path.exists(fs)):
            print(f'  {y}: missing temp/salt, skipped')
            continue
        with xr.open_dataset(ft, decode_times=False) as dt, \
             xr.open_dataset(fs, decode_times=False) as ds:
            T = dt['temp'].values
            S = ds['salt'].values
        if T.shape[0] != 12:
            print(f'  {y}: {T.shape[0]} time steps, not 12 -- skipped')
            continue
        this = np.empty((12, T.shape[1]))
        thisb = np.zeros((12, T.shape[1]))
        for m in range(12):
            mld, tob, _ = mld_woa_criterion(T[m].astype('f8'), S[m].astype('f8'), depth)
            this[m] = mld
            thisb[m] = tob.astype('f8')
        acc = this if acc is None else acc + this
        accb = thisb if accb is None else accb + thisb
        nyr += 1
        print(f'  {y}: done', flush=True)

        f2 = f'{root}/outdata/fesom/MLD2.fesom.{y}.nc'
        if os.path.exists(f2):
            with xr.open_dataset(f2, decode_times=False) as d2:
                nm = 'MLD2' if 'MLD2' in d2.data_vars else list(d2.data_vars)[-1]
                a2 = np.abs(d2[nm].values)               # FESOM writes it negative
            if a2.shape[0] == 12:
                acc2 = a2 if acc2 is None else acc2 + a2
    if nyr == 0:
        sys.exit('no usable model years found')
    return acc / nyr, (acc2 / nyr if acc2 is not None else None), accb / nyr, nyr


# -------------------------------------------------------------------------- WOA fields
def woa_monthly():
    """WOA18 M_an and M_dd, (12, 180, 360), plus the 1 degree axes."""
    mld = np.full((12, 180, 360), np.nan)
    dd = np.zeros((12, 180, 360))
    lat = lon = None
    for m in range(12):
        f = f'{WOA_DIR}/woa18_decav81B0_M02{m + 1:02d}_01.nc'
        if not os.path.exists(f):
            print(f'  WOA month {m + 1}: missing')
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            mld[m] = np.squeeze(d['M_an'].values)
            if 'M_dd' in d:
                dd[m] = np.nan_to_num(np.squeeze(d['M_dd'].values))
            lat = d['lat'].values
            lon = d['lon'].values
    return mld, dd, lat, lon


# ------------------------------------------------------------------------- regridding
def build_binner(mlon, mlat, wlat, wlon):
    """Node -> 1 degree cell index, and a nearest-node tree for cells with no node."""
    ilon = np.clip((mlon + 180.0).astype(int), 0, 359)
    ilat = np.clip((mlat + 90.0).astype(int), 0, 179)
    flat = ilat * 360 + ilon

    xyz = np.column_stack([
        np.cos(np.deg2rad(mlat)) * np.cos(np.deg2rad(mlon)),
        np.cos(np.deg2rad(mlat)) * np.sin(np.deg2rad(mlon)),
        np.sin(np.deg2rad(mlat))])
    tree = cKDTree(xyz)
    LO, LA = np.meshgrid(wlon, wlat)
    cxyz = np.column_stack([
        (np.cos(np.deg2rad(LA)) * np.cos(np.deg2rad(LO))).ravel(),
        (np.cos(np.deg2rad(LA)) * np.sin(np.deg2rad(LO))).ravel(),
        np.sin(np.deg2rad(LA)).ravel()])
    _, nearest = tree.query(cxyz, k=1)
    return flat, nearest


def to_grid(field, flat, area, nearest):
    """Area-weighted mean of nodes per 1 degree cell; nearest node where a cell has none."""
    ok = ~np.isnan(field)
    num = np.bincount(flat[ok], weights=(field[ok] * area[ok]), minlength=180 * 360)
    den = np.bincount(flat[ok], weights=area[ok], minlength=180 * 360)
    out = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
    empty = den == 0
    out[empty] = field[nearest[empty]]
    return out.reshape(180, 360)


# ----------------------------------------------------------------------- zonal, native
def zonal_native(field, lat, area, edges):
    """Band means in the grid's own geometry -- no regridding."""
    out = []
    for lo, hi in edges:
        s = (lat >= lo) & (lat < hi) & ~np.isnan(field)
        out.append(np.sum(field[s] * area[s]) / np.sum(area[s]) if s.any() else np.nan)
    return np.array(out)


def zonal_woa(field, wlat, edges):
    w = np.cos(np.deg2rad(wlat))[:, None] * np.ones((1, 360))
    out = []
    for lo, hi in edges:
        rows = (wlat >= lo) & (wlat < hi)
        f = field[rows]
        ww = w[rows]
        s = ~np.isnan(f)
        out.append(np.sum(f[s] * ww[s]) / np.sum(ww[s]) if s.any() else np.nan)
    return np.array(out)


# ------------------------------------------------------------------------------ plots
MLD_LEVELS = [0, 10, 20, 30, 40, 60, 80, 100, 150, 200, 300, 400, 600, 800, 1200, 2000]
BIAS_LEVELS = [-400, -300, -200, -150, -100, -60, -30, -10, 10, 30, 60, 100, 150, 200,
               300, 400]


def map_triptych(mod, obs, wlat, wlon, title, fname):
    fig, axes = plt.subplots(3, 1, figsize=(11, 13.5),
                             subplot_kw={'projection': ccrs.Robinson(central_longitude=200)})
    cmap = plt.get_cmap('viridis_r')
    norm = BoundaryNorm(MLD_LEVELS, cmap.N, extend='max')
    bcmap = plt.get_cmap('RdBu_r')
    bnorm = BoundaryNorm(BIAS_LEVELS, bcmap.N, extend='both')
    bias = mod - obs

    for ax, dat, lab, cm, nm, lv in (
            (axes[0], mod, f'AWI-ESM3 {ARM} ({YEARS[0]}-{YEARS[-1]}), WOA criterion',
             cmap, norm, MLD_LEVELS),
            (axes[1], obs, 'WOA18 M_an (1981-2010)', cmap, norm, MLD_LEVELS),
            (axes[2], bias, 'bias (model - obs); red = model mixes too deep',
             bcmap, bnorm, BIAS_LEVELS)):
        h = ax.pcolormesh(wlon, wlat, dat, transform=ccrs.PlateCarree(),
                          cmap=cm, norm=nm, shading='auto')
        ax.add_feature(cfeature.LAND, facecolor='0.85', zorder=2)
        ax.coastlines(linewidth=0.4, zorder=3)
        ax.set_global()
        ax.set_title(lab, fontsize=11)
        cb = fig.colorbar(h, ax=ax, orientation='vertical', pad=0.02, shrink=0.9,
                          ticks=lv[::2])
        cb.set_label('m', fontsize=9)
        cb.ax.tick_params(labelsize=8)
    fig.suptitle(title, fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    p = os.path.join(PLOTS, fname)
    fig.savefig(p, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {p}')


def seasonal_and_definition_plot(zm_woa_crit, zm_mld2, zm_obs, fname):
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharex=True)
    for i, (name, lo, hi) in enumerate(BANDS):
        ax = axes.ravel()[i]
        x = np.arange(1, 13)
        ax.plot(x, zm_obs[:, i], 'k-o', ms=3.5, lw=1.8, label='WOA18 M_an')
        ax.plot(x, zm_woa_crit[:, i], '-o', color='#c0392b', ms=3.5, lw=1.8,
                label=f'{ARM}, WOA criterion')
        if zm_mld2 is not None:
            ax.plot(x, zm_mld2[:, i], '--s', color='#2980b9', ms=3, lw=1.4,
                    label=f'{ARM}, FESOM MLD2 (surface ref)')
        ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.3)
        ax.invert_yaxis()
        ax.set_xticks(x)
        ax.set_xticklabels([m[0] for m in MONTHS], fontsize=8)
        if i % 4 == 0:
            ax.set_ylabel('MLD [m]  (deeper down)')
    ax = axes.ravel()[7]
    ax.axis('off')
    h, l = axes.ravel()[0].get_legend_handles_labels()
    ax.legend(h, l, loc='center', fontsize=10, frameon=False)
    ax.text(0.5, 0.18, 'band means computed on each grid\nin its own geometry '
                       '(no regridding)', ha='center', fontsize=9, color='0.35',
            transform=ax.transAxes)
    fig.suptitle(f'Mixed layer depth seasonal cycle by band: {ARM} '
                 f'({YEARS[0]}-{YEARS[-1]}) against WOA18 (1981-2010)', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(PLOTS, fname)
    fig.savefig(p, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {p}')


# ------------------------------------------------------------------------------- main
def main():
    os.makedirs(PLOTS, exist_ok=True)
    os.makedirs(os.path.join(REPO, 'data'), exist_ok=True)
    root = f'{R092}/{ARM}'
    if not os.path.isdir(root):
        cand = [p for p in glob.glob(f'{R092}/*{ARM}*') if os.path.isdir(p)]
        if not cand:
            sys.exit(f'no run directory for {ARM}')
        root = cand[0]

    print(f'ARM {ARM}  years {YEARS[0]}-{YEARS[-1]}  root {root}')
    mlon, mlat, area, depth = load_mesh()
    print(f'mesh ok: {mlon.size} nodes, {depth.size} levels, top mid-depths {depth[:4]}')

    print('model MLD (WOA criterion, per month per year):')
    mod, mld2, tob, nyr = model_monthly_mld(root, depth)
    print(f'  {nyr} years averaged')

    print('WOA18:')
    obs, dd, wlat, wlon = woa_monthly()

    flat, nearest = build_binner(mlon, mlat, wlat, wlon)

    mod_ann = np.nanmean(mod, axis=0)
    obs_ann = np.nanmean(obs, axis=0)
    mod_max = np.nanmax(mod, axis=0)
    obs_max = np.nanmax(obs, axis=0)


    mod_ann_g = to_grid(mod_ann, flat, area, nearest)
    mod_max_g = to_grid(mod_max, flat, area, nearest)

    # NOMAP=1 skips the cartopy maps (the system cartopy fails with a newer matplotlib)
    if os.environ.get('NOMAP', '0') != '1':
        map_triptych(mod_ann_g, obs_ann, wlat, wlon,
                     f'Annual-mean mixed layer depth  --  {ARM} vs WOA18, both on the '
                     f'0.125 kg/m$^3$ / 10 m criterion',
                     f'mld_annmean_{ARM}_vs_woa18.png')
        map_triptych(mod_max_g, obs_max, wlat, wlon,
                     f'Deepest month of the climatology (winter mixing)  --  {ARM} vs WOA18, '
                     f'same criterion',
                     f'mld_wintermax_{ARM}_vs_woa18.png')

    edges = [(lo, hi) for _, lo, hi in BANDS]
    zm_mod = np.array([zonal_native(mod[m], mlat, area, edges) for m in range(12)])
    zm_obs = np.array([zonal_woa(obs[m], wlat, edges) for m in range(12)])
    zm_m2 = (np.array([zonal_native(mld2[m], mlat, area, edges) for m in range(12)])
             if mld2 is not None else None)
    seasonal_and_definition_plot(zm_mod, zm_m2, zm_obs, f'mld_seasonal_bands_{ARM}.png')

    print('\n' + '=' * 96)
    print(f'MIXED LAYER DEPTH [m], {ARM} {YEARS[0]}-{YEARS[-1]} vs WOA18 1981-2010')
    print("both columns on WOA's criterion: sigma_theta(z)-sigma_theta(10m) >= 0.125 kg/m3")
    print('=' * 96)
    print(f'{"band":<24}{"ann model":>10}{"ann obs":>10}{"ann bias":>10}'
          f'{"max model":>11}{"max obs":>10}{"max bias":>10}')
    zmax_mod = zonal_native(mod_max, mlat, area, edges)
    zmax_obs = zonal_woa(obs_max, wlat, edges)
    zann_mod = zonal_native(mod_ann, mlat, area, edges)
    zann_obs = zonal_woa(obs_ann, wlat, edges)
    for i, (name, _, _) in enumerate(BANDS):
        print(f'{name:<24}{zann_mod[i]:>10.1f}{zann_obs[i]:>10.1f}'
              f'{zann_mod[i] - zann_obs[i]:>+10.1f}'
              f'{zmax_mod[i]:>11.1f}{zmax_obs[i]:>10.1f}'
              f'{zmax_mod[i] - zmax_obs[i]:>+10.1f}')

    if zm_m2 is not None:
        print('\n' + '-' * 96)
        print('DEFINITION SENSITIVITY: same model fields, two definitions [m, annual mean]')
        print(f'{"band":<24}{"WOA criterion":>16}{"FESOM MLD2":>14}{"difference":>13}')
        z2 = zm_m2.mean(axis=0)
        z1 = zm_mod.mean(axis=0)
        for i, (name, _, _) in enumerate(BANDS):
            print(f'{name:<24}{z1[i]:>16.1f}{z2[i]:>14.1f}{z1[i] - z2[i]:>+13.1f}')
        print('FESOM MLD2 references the surface level (2.5 m), WOA references 10 m.')

    print('\n' + '-' * 96)
    print('WOA OBSERVATIONAL COVERAGE: % of analysed ocean cells holding >=1 profile')
    print(f'{"band":<24}{"annual":>9}{"Mar":>7}{"Sep":>7}{"min month":>11}')
    for i, (name, lo, hi) in enumerate(BANDS):
        rows = (wlat >= lo) & (wlat < hi)
        oc = np.isfinite(obs[0][rows])
        n = max(oc.sum(), 1)
        cov = [(dd[m][rows][oc] >= 1).sum() / n * 100 for m in range(12)]
        print(f'{name:<24}{np.mean(cov):>9.1f}{cov[2]:>7.1f}{cov[8]:>7.1f}{min(cov):>11.1f}')

    fb = zonal_native(tob.mean(axis=0), mlat, area, edges)
    print('\nfraction of month-nodes mixed all the way to the bottom (model):')
    for i, (name, _, _) in enumerate(BANDS):
        print(f'  {name:<24}{fb[i]:>8.4f}')

    np.savez(os.path.join(REPO, 'data', f'mld_{ARM}_vs_woa18.npz'),
             mod_ann_grid=mod_ann_g, mod_max_grid=mod_max_g,
             obs_monthly=obs, obs_dd=dd, wlat=wlat, wlon=wlon,
             zm_mod=zm_mod, zm_obs=zm_obs,
             zm_mld2=(zm_m2 if zm_m2 is not None else np.array(np.nan)),
             years=np.array(YEARS))
    print(f'\nsaved data/mld_{ARM}_vs_woa18.npz')


if __name__ == '__main__':
    main()
