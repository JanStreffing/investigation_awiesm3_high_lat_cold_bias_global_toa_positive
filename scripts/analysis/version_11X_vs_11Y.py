"""11Y against 11X: do FESOM main's ~273 commits move this configuration?

THE TEST.  One variable, the FESOM library.  11X stages libfesom.so built 2026-08-26 from
awiesm3-implicit-ice-surftemp; 11Y stages the 09-08 build from main.  Everything else is
byte-identical: same corrected CORE3 mesh, same core3_linfs restart, same LPJ-GUESS binary
and state, same 1792 OASIS set, same namelists.  Both ran their four legs on ONE library
throughout, so neither arm is internally mixed.

WHY IT COULD MATTER.  main carries EVP mass-regularisation (#951), EVP frozen stress (#950),
a KPP Stokes divide-by-zero fix (#948) and a mass-matrix tolerance fix (#1043) -- ice
dynamics and vertical mixing, the two things the Weddell polynya depends on.

A NULL IS THE USEFUL ANSWER HERE, and it is the expected one.  If nothing exceeds the
detection threshold, the FESOM version can be updated without re-baselining the campaign,
and 11X-era results stay comparable with anything built on main.  If something does move,
every cross-version comparison in the campaign needs revisiting.

MESH GUARD.  Both arms are the corrected core3 (220509 nodes) -- unlike 11X vs 11W, this is
NOT a mesh pair.  The guard is kept and asserted anyway, because loading the wrong mesh
fails silently through numpy broadcasting.

TRAPS OBSERVED (campaign protocol, unchanged).
    IFS accumulated fluxes are J/m2 over the output step; divide by 3600.
    Southern Ocean band means MUST be ocean-masked; a 90-60S average folds in the ice sheet
    and halved a polar signal once already.
    Sea-ice AREA is not EXTENT.  Both reported, separately labelled.
    Thresholds are the control's own interannual scatter, 1.96*sd*sqrt(2/n), paired form.
    Band means are computed on each grid in its own geometry - no regridding.

CAVEAT CARRIED.  Both arms ran with the ice-skin defects found on 2026-09-08.  They share
them, so a DIFFERENCE between the arms is readable; absolute values are not clean.  Neither
arm has the MLE scheme: 12A-12D are not comparable to these.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
WOA = '/work/ab0246/a270092/obs/WOA/woa18_decav81B0_M0216_01.nc'
LSMF = ('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
ACC = 3600.0

ARMS = [('11X fesom 08-26', f'{R092}/11X', '/work/ab0246/a270092/input/fesom2/core3', 220509),
        ('11Y fesom main',  f'{R092}/11Y', '/work/ab0246/a270092/input/fesom2/core3', 220509)]

STATE = list(range(1382, 1390))     # last 8 common years, as mld_vs_woa18.py used
DRIFT = list(range(1350, 1390))     # common span
FALL = [9, 10, 11]                  # 0-based Oct, Nov, Dec, to match WOA18 fall
BANDS = [('90-60S  Antarctic/SO', -90, -60), ('60-45S  subantarctic', -60, -45),
         ('45-30S', -45, -30), ('30S-30N  tropics', -30, 30),
         ('45-60N  subpolar NA', 45, 60), ('60-90N  Arctic/Nordic', 60, 90)]
BANDS = BANDS[::-1]          # report north to south
SIB = (55.0, 75.0, 60.0, 180.0)
DJF, JJA = [11, 0, 1], [5, 6, 7]    # 0-based month indices


def mesh_of(path, nexp):
    """Node latitude and node area, from whichever file this mesh ships."""
    diag = f'{path}/fesom.mesh.diag.nc'
    if os.path.exists(diag):
        with xr.open_dataset(diag, decode_times=False) as m:
            lat = m['lat'].values
            area = m['nod_area'].values[0, :]
    else:
        with xr.open_dataset(f'{path}/mesh.nc', decode_times=False) as m:
            lat = m['lat'].values
            area = m['cell_area'].values
    if np.abs(lat).max() < 4:
        lat = np.rad2deg(lat)
    assert lat.size == nexp, f'MESH GUARD: {path} has {lat.size} nodes, expected {nexp}'
    return lat, area


def fesom_years(root, var, years, months=None):
    """Yield (year, field) with the month selection applied, native grid."""
    for y in years:
        f = glob.glob(f'{root}/outdata/fesom/{var}.fesom.{y}.nc')
        if not f:
            continue
        with xr.open_dataset(f[0], decode_times=False) as d:
            nm = var if var in d.data_vars else list(d.data_vars)[-1]
            a = d[nm].values
        yield y, (a[months] if months is not None else a)


def oifs_annual(root, var, years):
    """Annual means of a remapped monthly OIFS field, plus lat/lon."""
    out, lat, lon = [], None, None
    for y in years:
        f = glob.glob(f'{root}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc')
        if not f:
            continue
        with xr.open_dataset(f[0], decode_times=False) as d:
            a = d[var].values
            lat, lon = d['lat'].values, d['lon'].values
        out.append((y, a))
    return out, lat, lon


def band(vals, lat, area, la0, la1):
    m = (lat >= la0) & (lat < la1) & np.isfinite(vals)
    return float(np.sum(vals[m] * area[m]) / np.sum(area[m])) if m.any() else np.nan


def thr(series):
    """Paired detection threshold from interannual scatter: 1.96*sd*sqrt(2/n)."""
    s = np.asarray([v for v in series if np.isfinite(v)])
    return 1.96 * s.std(ddof=1) * np.sqrt(2.0 / s.size) if s.size > 2 else np.nan


print(__doc__.split('\n')[0])
print(f'state window {STATE[0]}-{STATE[-1]}, drift window {DRIFT[0]}-{DRIFT[-1]}\n')

# ---------------------------------------------------------------- 1. the falsifiers
print('=' * 78)
print('1. THE PRE-REGISTERED TEST: Weddell convection')
print('=' * 78)
with xr.open_dataset(WOA, decode_times=False) as d:
    wmld = np.squeeze(d['M_an'].values); wlat = d['lat'].values
ws = (wlat >= -90) & (wlat < -60)
sub = wmld[ws, :]
wgt = np.broadcast_to(np.cos(np.deg2rad(wlat[ws]))[:, None], sub.shape)
ok = np.isfinite(sub)
woa_polar = float(np.average(sub[ok], weights=wgt[ok]))
print(f'  WOA18 Oct-Dec 90-60S mixed layer: {woa_polar:.1f} m   (the target)\n')

print(f'  {"arm":15s} {"MLD2 90-60S":>12s} {"vs WOA18":>10s} {">1000m area":>12s} {"deepest":>9s}')
print('  ' + '-' * 62)
mld_store = {}
for label, root, mpath, nexp in ARMS:
    lat, area = mesh_of(mpath, nexp)
    acc, n = None, 0
    for y, a in fesom_years(root, 'MLD2', STATE, FALL):
        assert a.shape[-1] == nexp, f'MESH GUARD: {label} MLD2 has {a.shape[-1]} nodes'
        v = np.abs(a).mean(axis=0)
        acc = v if acc is None else acc + v; n += 1
    if not n:
        print(f'  {label:15s} no MLD2'); continue
    m2 = acc / n
    pol = band(m2, lat, area, -90, -60)
    sel = (lat >= -90) & (lat < -60)
    frac = float(np.sum(area[sel][m2[sel] > 1000.0]) / np.sum(area[sel])) * 100.0
    mld_store[label] = (m2, lat, area)
    print(f'  {label:15s} {pol:12.1f} {pol - woa_polar:+10.1f} {frac:11.1f}% {np.nanmax(m2[sel]):9.0f}')

print(f'\n  {"band":24s} {"WOA18":>9s} ' + ' '.join(f'{a[0]:>15s}' for a in ARMS))
print('  ' + '-' * 72)
for name, la0, la1 in BANDS:
    bs = (wlat >= la0) & (wlat < la1)
    s = wmld[bs, :]
    w = np.broadcast_to(np.cos(np.deg2rad(wlat[bs]))[:, None], s.shape)
    o = np.isfinite(s)
    wo = float(np.average(s[o], weights=w[o])) if o.any() else np.nan
    row = f'  {name:24s} {wo:9.1f} '
    for label, *_ in ARMS:
        if label in mld_store:
            m2, lat, area = mld_store[label]
            row += f'{band(m2, lat, area, la0, la1):15.1f} '
    print(row)

# ------------------------------------------------- 1b. is that difference resolvable?
print('\n' + '=' * 78)
print('1b. SIGNIFICANCE of the MLD change (per-year, paired threshold)')
print('=' * 78)
peryear = {}
for label, root, mpath, nexp in ARMS:
    lat, area = mesh_of(mpath, nexp)
    sel = (lat >= -90) & (lat < -60)
    vals, fracs = [], []
    for y, a in fesom_years(root, 'MLD2', DRIFT, FALL):
        v = np.abs(a).mean(axis=0)
        vals.append((y, band(v, lat, area, -90, -60)))
        fracs.append((y, float(np.sum(area[sel][v[sel] > 1000.0]) / np.sum(area[sel])) * 100.0))
    peryear[label] = (vals, fracs)
    s = [v for _, v in vals if np.isfinite(v)]
    print(f'  {label:15s} n={len(s):2d} yr   mean {np.mean(s):7.1f} m   sd {np.std(s, ddof=1):6.1f}   '
          f'range {min(s):.0f}-{max(s):.0f}')

(la, lb) = [a[0] for a in ARMS]
sa = np.array([v for _, v in peryear[la][0]]); sb = np.array([v for _, v in peryear[lb][0]])
n = min(sa.size, sb.size)
t_mld = thr(sa[-8:])                      # control scatter over the state window
d_mld = np.mean(sb[-8:]) - np.mean(sa[-8:])
fa = np.array([v for _, v in peryear[la][1]]); fb = np.array([v for _, v in peryear[lb][1]])
t_fr = thr(fa[-8:]); d_fr = np.mean(fb[-8:]) - np.mean(fa[-8:])
print(f'\n  MLD2 90-60S   11X - 11W = {d_mld:+7.1f} m    threshold +-{t_mld:.1f}    '
      f'{"RESOLVED" if abs(d_mld) > t_mld else "within noise"}')
print(f'  >1000m area   11X - 11W = {d_fr:+7.1f} pp   threshold +-{t_fr:.1f}    '
      f'{"RESOLVED" if abs(d_fr) > t_fr else "within noise"}')
gap_before = np.mean(sa[-8:]) - woa_polar
print(f'\n  gap to WOA18: {gap_before:+.1f} m before, {np.mean(sb[-8:]) - woa_polar:+.1f} m after '
      f'-> {100 * (-d_mld) / gap_before:.1f}% of the bias removed')

# --------------------------- 1c. is the last-8 window representative, or was it quiet?
print('\n' + '=' * 78)
print('1c. THRESHOLD SENSITIVITY: which scatter do we trust?')
print('=' * 78)
sd8, sd40 = sa[-8:].std(ddof=1), sa.std(ddof=1)
t8 = 1.96 * sd8 * np.sqrt(2 / 8); t40 = 1.96 * sd40 * np.sqrt(2 / 8)
print(f'  control (11W) sd over last 8 yr : {sd8:6.1f} m  -> paired threshold +-{t8:6.1f} m')
print(f'  control (11W) sd over all 40 yr : {sd40:6.1f} m  -> paired threshold +-{t40:6.1f} m')
print(f'  measured difference             : {d_mld:+6.1f} m')
print(f'  verdict on last-8 scatter : {"RESOLVED" if abs(d_mld) > t8 else "within noise"}')
print(f'  verdict on full-run scatter: {"RESOLVED" if abs(d_mld) > t40 else "WITHIN NOISE"}')

print('\n  decadal means, 90-60S Oct-Dec MLD2 [m]:')
print(f'  {"decade":12s} ' + ' '.join(f'{a[0]:>15s}' for a in ARMS) + '     diff')
for d0 in (1350, 1360, 1370, 1380):
    row = f'  {d0}-{d0+9:4d} '
    vals = []
    for label, *_ in ARMS:
        v = [x for y, x in peryear[label][0] if d0 <= y < d0 + 10 and np.isfinite(x)]
        vals.append(np.mean(v) if v else np.nan)
        row += f'{vals[-1]:15.1f} '
    row += f'  {vals[1] - vals[0]:+7.1f}'
    print(row)

# full-run trend: is either arm drifting?
for label, *_ in ARMS:
    ys = np.array([y for y, v in peryear[label][0] if np.isfinite(v)], dtype=float)
    vs = np.array([v for _, v in peryear[label][0] if np.isfinite(v)])
    sl = np.polyfit(ys, vs, 1)[0]
    print(f'  {label:15s} trend {sl * 10:+7.1f} m/decade')

print('\n  >1000 m convecting area fraction [%], same treatment:')
sd8f, sd40f = fa[-8:].std(ddof=1), fa.std(ddof=1)
print(f'    control sd  last-8 {sd8f:5.2f} -> +-{1.96*sd8f*np.sqrt(2/8):5.2f} pp ; '
      f'all-40 {sd40f:5.2f} -> +-{1.96*sd40f*np.sqrt(2/8):5.2f} pp ; measured {d_fr:+.2f} pp')
row = '    decadal diff (11X-11W): '
for d0 in (1350, 1360, 1370, 1380):
    va = np.mean([x for y, x in peryear[la][1] if d0 <= y < d0 + 10])
    vb = np.mean([x for y, x in peryear[lb][1] if d0 <= y < d0 + 10])
    row += f'{d0}s {vb - va:+6.2f}  '
print(row)

# ------------------------------------------- 2. third falsifier: the polar halocline
print('\n' + '=' * 78)
print('2. UPPER-50 m SALINITY at 90-60S (the halocline the convection destroys)')
print('=' * 78)
print('  The PHC3 reference on disk is core3_beta only (211567 nodes), so a BIAS can be')
print('  quoted for 11W and not for 11X.  Absolutes are comparable: each is a band mean')
print('  in its own geometry, no regridding.')
PHC3 = '/work/ab0246/a270092/postprocessing/climatologies/CORE3/salt.fesom.1958.nc'
with xr.open_dataset(PHC3, decode_times=False) as r:
    ref_s = np.squeeze(r['salt'].values)              # (nod2, nz1) on core3_beta
with xr.open_dataset('/work/ab0246/a270092/input/fesom2/core3_beta/fesom.mesh.diag.nc',
                     decode_times=False) as m:
    zmid = m['nz1'].values
ktop = np.where(zmid <= 50.0)[0]
print(f'\n  upper-50 m = layers {ktop[0]}..{ktop[-1]} '
      f'(mid-depths {zmid[ktop][0]:.1f}-{zmid[ktop][-1]:.1f} m)')
print(f'\n  {"arm":15s} {"S 90-60S":>10s} {"vs PHC3":>10s}')
print('  ' + '-' * 40)
salt_abs = {}
for label, root, mpath, nexp in ARMS:
    lat, area = mesh_of(mpath, nexp)
    acc, n = None, 0
    for y, a in fesom_years(root, 'salt', STATE):        # (time, nod2, nz)
        assert a.shape[1] == nexp, f'MESH GUARD: {label} salt has {a.shape[1]} nodes'
        v = a.mean(axis=0)[:, ktop].mean(axis=1)
        acc = v if acc is None else acc + v; n += 1
    if not n:
        print(f'  {label:15s} no salt'); continue
    val = band(acc / n, lat, area, -90, -60)
    salt_abs[label] = val
    if nexp == 211567:
        rb = band(ref_s[:, ktop].mean(axis=1), lat, area, -90, -60)
        print(f'  {label:15s} {val:10.3f} {val - rb:+10.3f}   (PHC3 {rb:.3f})')
    else:
        print(f'  {label:15s} {val:10.3f} {"n/a":>10s}   (no core3 PHC3 on disk)')
if len(salt_abs) == 2:
    print(f'\n  11X - 11W = {salt_abs[lb] - salt_abs[la]:+.3f} psu '
          f'(negative = fresher = halocline restored)')

# ---------------------------------------------- 3. the standard coupled metrics
print('\n' + '=' * 78)
print('3. STANDARD COUPLED METRICS  (campaign protocol)')
print('=' * 78)
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values)[0]
ocean = lsm < 0.5

def toa_series(root, years):
    """Net TOA [W/m2] per year, and SW CRE, on the remapped OIFS grid."""
    net, cre, t2m, ci_sh = [], [], [], []
    for y in years:
        got = {}
        for v in ('tsr', 'ttr', 'tsrc', '2t', 'ci'):
            f = glob.glob(f'{root}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc')
            if f:
                with xr.open_dataset(f[0], decode_times=False) as d:
                    got[v] = d[v].values
                    lat_, lon_ = d['lat'].values, d['lon'].values
        if 'tsr' not in got or 'ttr' not in got:
            continue
        w = np.cos(np.deg2rad(lat_))[:, None] * np.ones((1, lat_.size and got['tsr'].shape[-1]))
        n = ((got['tsr'] + got['ttr']) / ACC).mean(axis=0)
        net.append((y, float(np.average(n, weights=w))))
        if 'tsrc' in got:
            c = ((got['tsr'] - got['tsrc']) / ACC).mean(axis=0)
            cre.append((y, c, lat_))
        if '2t' in got:
            t2m.append((y, got['2t'], lat_, lon_))
        if 'ci' in got:
            ci_sh.append((y, got['ci'], lat_))
    return net, cre, t2m, ci_sh

print(f'  {"arm":15s} {"net TOA":>9s} {"drift":>12s} {"SibJJA T2m":>11s} {"SibDJF T2m":>11s}')
print('  ' + '-' * 62)
store = {}
for label, root, mpath, nexp in ARMS:
    net, cre, t2m, ci_sh = toa_series(root, DRIFT)
    ys = np.array([y for y, _ in net], float); vs = np.array([v for _, v in net])
    sl = np.polyfit(ys, vs, 1)[0] * 10 if ys.size > 2 else np.nan
    la0, la1, lo0, lo1 = SIB
    sj = sd_ = np.nan
    if t2m:
        lat_, lon_ = t2m[0][2], t2m[0][3]
        ms = (lat_ >= la0) & (lat_ < la1); mo = (lon_ >= lo0) & (lon_ < lo1)
        w = np.cos(np.deg2rad(lat_[ms]))[:, None] * np.ones((1, mo.sum()))
        jj = [float(np.average(a[JJA][:, ms][:, :, mo].mean(axis=0), weights=w))
              for y, a, *_ in t2m if y >= STATE[0]]
        dj = [float(np.average(a[DJF][:, ms][:, :, mo].mean(axis=0), weights=w))
              for y, a, *_ in t2m if y >= STATE[0]]
        sj, sd_ = np.mean(jj), np.mean(dj)
    store[label] = dict(net=vs, ys=ys, cre=cre, ci=ci_sh, jj=sj, dj=sd_)
    print(f'  {label:15s} {vs[-8:].mean():9.3f} {sl:+9.3f}/dec {sj - 273.15:11.2f} {sd_ - 273.15:11.2f}')

# per-year series for thresholds
print('\n  significance, thresholds from the CONTROL (11W) interannual scatter:')
print(f'  {"metric":22s} {"11W":>9s} {"11X":>9s} {"diff":>9s} {"thr(8y)":>9s} {"thr(40y)":>9s}  verdict')
print('  ' + '-' * 84)

def sib_series(root, months):
    out = []
    for y in STATE:
        f = glob.glob(f'{root}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc')
        if not f:
            continue
        with xr.open_dataset(f[0], decode_times=False) as d:
            a = d['2t'].values; lat_ = d['lat'].values; lon_ = d['lon'].values
        ms = (lat_ >= SIB[0]) & (lat_ < SIB[1]); mo = (lon_ >= SIB[2]) & (lon_ < SIB[3])
        w = np.cos(np.deg2rad(lat_[ms]))[:, None] * np.ones((1, mo.sum()))
        out.append(float(np.average(a[months][:, ms][:, :, mo].mean(axis=0), weights=w)))
    return np.array(out)

rows = []
rows.append(('net TOA [W/m2]', store[la]['net'][-8:], store[lb]['net'][-8:], store[la]['net']))
rows.append(('Siberian JJA T2m [K]', sib_series(ARMS[0][1], JJA), sib_series(ARMS[1][1], JJA), None))
rows.append(('Siberian DJF T2m [K]', sib_series(ARMS[0][1], DJF), sib_series(ARMS[1][1], DJF), None))
for name, a8, b8, full in rows:
    d = b8.mean() - a8.mean()
    t8 = 1.96 * a8.std(ddof=1) * np.sqrt(2 / a8.size)
    t40 = 1.96 * full.std(ddof=1) * np.sqrt(2 / a8.size) if full is not None else np.nan
    tref = t40 if np.isfinite(t40) else t8
    v = 'RESOLVED' if abs(d) > tref else 'within noise'
    print(f'  {name:22s} {a8.mean():9.3f} {b8.mean():9.3f} {d:+9.3f} {t8:9.3f} '
          f'{t40:9.3f}  {v}')

# sea ice: AREA and EXTENT are different quantities
print('\n  sea ice, September SH and March NH (area = sum a*A; extent = sum A where a>0.15):')
print(f'  {"arm":15s} {"SH Sep area":>12s} {"SH Sep ext":>12s} {"NH Mar area":>12s} {"NH Mar ext":>12s}')
print('  ' + '-' * 68)
for label, root, mpath, nexp in ARMS:
    lat, area = mesh_of(mpath, nexp)
    sh, nh = lat < 0, lat > 0
    acc = {k: [] for k in ('sa', 'se', 'na', 'ne')}
    for y, a in fesom_years(root, 'a_ice', STATE):
        assert a.shape[-1] == nexp, f'MESH GUARD: {label} a_ice {a.shape[-1]} nodes'
        sep, mar = a[8], a[2]
        acc['sa'].append(np.sum(sep[sh] * area[sh]) / 1e12)
        acc['se'].append(np.sum(area[sh][sep[sh] > 0.15]) / 1e12)
        acc['na'].append(np.sum(mar[nh] * area[nh]) / 1e12)
        acc['ne'].append(np.sum(area[nh][mar[nh] > 0.15]) / 1e12)
    print(f'  {label:15s} {np.mean(acc["sa"]):12.2f} {np.mean(acc["se"]):12.2f} '
          f'{np.mean(acc["na"]):12.2f} {np.mean(acc["ne"]):12.2f}   [1e6 km2]')

# ---- Siberian T2m on the FULL run, because the 8-year window misled on MLD
print('\n  Siberian T2m over the full 40 years (the 8-year window understated MLD scatter):')
def sib_full(root, months):
    out = []
    for y in DRIFT:
        f = glob.glob(f'{root}/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc')
        if not f:
            continue
        with xr.open_dataset(f[0], decode_times=False) as d:
            a = d['2t'].values; lat_ = d['lat'].values; lon_ = d['lon'].values
        ms = (lat_ >= SIB[0]) & (lat_ < SIB[1]); mo = (lon_ >= SIB[2]) & (lon_ < SIB[3])
        w = np.cos(np.deg2rad(lat_[ms]))[:, None] * np.ones((1, mo.sum()))
        out.append((y, float(np.average(a[months][:, ms][:, :, mo].mean(axis=0), weights=w))))
    return out
for nm, mo in (('JJA', JJA), ('DJF', DJF)):
    A = sib_full(ARMS[0][1], mo); B = sib_full(ARMS[1][1], mo)
    va = np.array([v for _, v in A]); vb = np.array([v for _, v in B])
    d8 = vb[-8:].mean() - va[-8:].mean()
    t40 = 1.96 * va.std(ddof=1) * np.sqrt(2 / 8)
    print(f'    Siberian {nm}: 11W sd(40y) {va.std(ddof=1):.3f} K -> thr +-{t40:.3f} ; '
          f'diff(last 8) {d8:+.3f} -> {"RESOLVED" if abs(d8) > t40 else "WITHIN NOISE"}')
    row = '      decadal diff: '
    for d0 in (1350, 1360, 1370, 1380):
        xa = np.mean([v for y, v in A if d0 <= y < d0 + 10])
        xb = np.mean([v for y, v in B if d0 <= y < d0 + 10])
        row += f'{d0}s {xb - xa:+6.2f}  '
    print(row)

print('\n  sea ice AREA, recomputed with nansum (a_ice carries fill values):')
for label, root, mpath, nexp in ARMS:
    lat, area = mesh_of(mpath, nexp)
    sh, nh = lat < 0, lat > 0
    sa_, na_ = [], []
    for y, a in fesom_years(root, 'a_ice', STATE):
        sep, mar = np.nan_to_num(a[8]), np.nan_to_num(a[2])
        sep = np.where(np.abs(sep) > 1.5, 0.0, sep); mar = np.where(np.abs(mar) > 1.5, 0.0, mar)
        sa_.append(np.sum(sep[sh] * area[sh]) / 1e12)
        na_.append(np.sum(mar[nh] * area[nh]) / 1e12)
    print(f'    {label:15s} SH Sep area {np.mean(sa_):6.2f}   NH Mar area {np.mean(na_):6.2f}  [1e6 km2]')

# ---- TOA drift is the campaign's energy target; quote its uncertainty
print('\n  net TOA drift with slope uncertainty (the energy target):')
for label, *_ in ARMS:
    ys, vs = store[label]['ys'], store[label]['net']
    n = ys.size
    b, a = np.polyfit(ys, vs, 1)
    resid = vs - (b * ys + a)
    se = np.sqrt((resid ** 2).sum() / (n - 2) / ((ys - ys.mean()) ** 2).sum())
    print(f'    {label:15s} {b*10:+6.3f} +- {1.96*se*10:.3f} W/m2/decade   '
          f'(level last-8 {vs[-8:].mean():+.3f})')
d_slope = (np.polyfit(store[lb]["ys"], store[lb]["net"], 1)[0]
           - np.polyfit(store[la]["ys"], store[la]["net"], 1)[0]) * 10
print(f'    difference in drift: {d_slope:+.3f} W/m2/decade')
