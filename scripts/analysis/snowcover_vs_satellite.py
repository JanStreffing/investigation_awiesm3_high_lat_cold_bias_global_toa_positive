"""Model snow cover vs SATELLITE, not vs another land model.

Everything in the round-13/14 snow story rested on ERA5, whose land surface is the
same HTESSEL family as ours -- so comparing our snow fields to ERA5's was partly
comparing a scheme to a sibling of itself. CERES settled the ALBEDO half
independently (its surface albedo is a radiative retrieval, and it confirms the
model's June excess at +0.046). This script settles the other half: snow COVER and
melt TIMING, against an actual satellite record.

Reference: Rutgers NH 24 km Weekly Snow Cover Extent CDR (NSIDC G10035), 1980-2024.
Visible-band NOAA snow charts through May 1999, NOAA/NIC IMS thereafter. Downloaded
to obs/snowcover/.

WHAT IS AND IS NOT COMPARABLE. Rutgers SCE is a BINARY classification per 24 km cell;
the box mean here is the area-weighted fraction of land cells classed snow-covered.
HTESSEL's ZCVS is a SUB-GRID fraction within a ~100 km cell. The two agree on
"how much of the region is snow" but differ in how a partially covered cell is
counted, which matters most mid-melt. The MELT TIMING comparison -- when does cover
fall through 50 %, and how fast -- is robust to that; the absolute mid-melt fraction
is less so, and is flagged rather than over-read.

Note also that ERA5 assimilates IMS snow cover, so ERA5 agreeing with Rutgers is
expected and is not independent corroboration. What it does establish is that
ERA5's extent is anchored to this observation, which is why extent was treated as
its trustworthy field all along.
"""
import numpy as np, xarray as xr, os, sys, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF, Y0, Y1

BOX = ((55, 75), (60, 180))
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
SAT = '/work/ab0246/a270092/obs/snowcover/rutgers_nh_24km_weekly_sce.nc'
E5 = '/work/ab0246/a270092/obs/era5/snow'
PD = list(range(1990, 2015))
PI = list(range(Y0, Y1 + 1))
RQSNCR = 1.0 / 10.0

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


# ------------------------------------------------------------------ satellite
if not os.path.exists(SAT):
    sys.exit(f'missing {SAT}')
d = xr.open_dataset(SAT)
sla = d['latitude'].values
slo = d['longitude'].values
land = d['land'].values
area = d['area'].values
t = d['time'].values                     # seconds since 1970 (week end date)

l180 = ((slo + 180) % 360) - 180
sel = ((sla >= BOX[0][0]) & (sla <= BOX[0][1]) &
       (l180 >= BOX[1][0]) & (l180 <= BOX[1][1]) & (land == 1))
print(f'Rutgers grid: {sel.sum()} land cells inside the Siberian box '
      f'({area[sel].sum()/1e6:.2f} million km2)')

tt = np.array(t, dtype='datetime64[s]')
yrs = tt.astype('datetime64[Y]').astype(int) + 1970
mos = tt.astype('datetime64[M]').astype(int) % 12 + 1
keep = (yrs >= PD[0]) & (yrs <= PD[-1])
w = area[sel]

sce = np.full(12, np.nan)
for m in range(12):
    idx = np.where(keep & (mos == m + 1))[0]
    if not len(idx):
        continue
    vals = []
    for i in idx:
        s = d['snow_cover_extent'].values[i][sel]
        ok = np.isfinite(s) & (s >= 0)
        vals.append(np.average(s[ok], weights=w[ok]) if ok.any() else np.nan)
    sce[m] = np.nanmean(vals)
d.close()


# ------------------------------------------------------------------ model / ERA5
def bm(a, lat, lon):
    ys = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
    ll = ((lon + 180) % 360) - 180
    xs = (ll >= BOX[1][0]) & (ll <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub = a[ii]
    mk = np.isfinite(sub) & (lsm[ii] > 0.5)
    ww = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape)
    return np.average(sub[mk], weights=ww[mk]) if mk.any() else np.nan


def clim(run, var, years):
    acc, n, lat, lon = None, 0, None, None
    for y in years:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        ds = xr.open_dataset(f)
        a = ds[var].values
        lat, lon = ds[var].lat.values, ds[var].lon.values
        ds.close()
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, lat, lon) if n else (None, None, None)


# Which snow-cover scheme each run actually used. Reconstructing f_snow with the
# as-released formula for a run that used Niu & Yang would report a cover the model
# never saw -- so the scheme is keyed by run, never assumed.
#   ('rqsncr', rq)          ZCVS = min(1, d_cm*rq)          as released
#   ('niuyang', z0, m, rn)  SCF  = tanh(d/(2.5*z0*(rho/rn)**m))
SCHEME = {
    'amip_H1_snowcr30':    ('rqsncr', 1.0 / 30.0),
    'amip_H2_G1_snowcr30': ('rqsncr', 1.0 / 30.0),
    'amip_I1_scf':         ('niuyang', 0.016, 1.6, 100.0),
    'amip_I2_scf_only':    ('niuyang', 0.016, 1.6, 100.0),
    'amip_I3_scf_sdor':    ('niuyang', 0.014, 1.6, 100.0),   # + CSD*SDOR, see note below
}


def model_fsnow(run, years, rq=RQSNCR):
    sd, lat, lon = clim(run, 'sd', years)
    rsn, _, _ = clim(run, 'rsn', years)
    if sd is None:
        return None
    sch = SCHEME.get(run, ('rqsncr', rq))
    out = []
    for m in range(12):
        d_m = sd[m] * 1000.0 / np.maximum(rsn[m], 1.0)          # depth [m]
        if sch[0] == 'rqsncr':
            fs = np.clip(d_m * 100.0 * sch[1], 0, 1)
        else:
            _, z0, mm, rn = sch
            scl = 2.5 * z0 * np.power(np.maximum(rsn[m], 50.0) / rn, mm)
            fs = np.clip(np.tanh(d_m / np.maximum(scl, 1e-9)), 0, 1)
        out.append(bm(fs, lat, lon))
    return np.array(out)


# ERA5, using the identical HTESSEL formula on ERA5's snow mass
dd = xr.open_dataset(f'{E5}/era5_141_clim_1990-2014.nc'); sd5 = dd['var141'].values
la, lo = dd['lat'].values, dd['lon'].values; dd.close()
dd = xr.open_dataset(f'{E5}/era5_033_clim_1990-2014.nc'); rsn5 = dd['var33'].values; dd.close()
dd = xr.open_dataset(f'{E5}/era5_lsm.nc'); lsm5 = np.squeeze(dd['var172'].values); dd.close()


def bm5(a):
    ys = (la >= BOX[0][0]) & (la <= BOX[0][1])
    ll = ((lo + 180) % 360) - 180
    xs = (ll >= BOX[1][0]) & (ll <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub = a[ii]
    mk = np.isfinite(sub) & (lsm5[ii] > 0.5)
    ww = np.broadcast_to(np.cos(np.deg2rad(la[ys]))[:, None], sub.shape)
    return np.average(sub[mk], weights=ww[mk]) if mk.any() else np.nan


e5f = np.array([bm5(np.clip(100.0 * sd5[m] * 1000.0 / np.maximum(rsn5[m], 1.0) * RQSNCR, 0, 1))
                for m in range(12)])

runs = [('presentday', 'amip_presentday', PD), ('control', 'amip_pi_base', PI),
        ('G4 tundra', 'amip_G4_tundra', PI),
        ('I1 scf', 'amip_I1_scf', PI), ('I2 scf only', 'amip_I2_scf_only', PI)]
M = {}
for lab, r, years in runs:
    v = model_fsnow(r, years)
    if v is not None:
        M[lab] = v

print(f'\nSiberian land box, SNOW COVER FRACTION -- satellite vs models\n')
print(f'  {"":5s} {"Rutgers":>9s} {"ERA5":>8s} ' + ' '.join(f'{l:>12s}' for l in M))
for m in range(12):
    row = ' '.join(f'{M[l][m]:>12.3f}' for l in M)
    mk = '  <<<' if m in (3, 4, 5) else ''
    print(f'  {MON[m]:5s} {sce[m]:9.3f} {e5f[m]:8.3f} {row}{mk}')

print(f'\n  differences vs SATELLITE (Rutgers):')
print(f'  {"":5s} {"ERA5":>9s} ' + ' '.join(f'{l:>12s}' for l in M))
for m in range(12):
    row = ' '.join(f'{M[l][m]-sce[m]:>+12.3f}' for l in M)
    mk = '  <<<' if m in (3, 4, 5) else ''
    print(f'  {MON[m]:5s} {e5f[m]-sce[m]:>+9.3f}{row}{mk}')

# melt timing: linear interpolation of the month at which cover falls through 0.5
def cross(v, thr=0.5):
    for m in range(3, 8):
        if v[m] >= thr > v[m + 1]:
            return m + (v[m] - thr) / max(v[m] - v[m + 1], 1e-9)
    return np.nan


print(f'\n  MELT TIMING -- month index at which snow cover falls through 50 % '
      f'(3.0 = 1 Apr, 4.0 = 1 May):')
print(f'    Rutgers (satellite) {cross(sce):.2f}')
print(f'    ERA5                {cross(e5f):.2f}')
for l in M:
    print(f'    {l:19s} {cross(M[l]):.2f}')
print('\n  Caveat: Rutgers is a binary 24 km classification, HTESSEL ZCVS is a sub-grid')
print('  fraction, so mid-melt absolute values are not strictly comparable. The')
print('  CROSSING DATE is the robust comparison.')
