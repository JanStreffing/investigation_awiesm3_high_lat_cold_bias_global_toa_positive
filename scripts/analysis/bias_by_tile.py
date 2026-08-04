"""Seasonal T2m bias against ERA5, resolved by dominant surface type.

Two questions this answers.

1. IS THE COLD BIAS A BOREAL PHENOMENON, or does the model run cold over most
   land? The campaign has only ever scored one box (Siberia 55-75N, 60-180E), so
   it has never checked whether the levers it tunes there are fixing a local
   fault or a global one. A lever indexed by VEGETATION TYPE acts wherever that
   type occurs, so this matters for every surface lever, not just the current one.

2. DOES A TYPE'S NON-BOREAL AREA SHARE THE BIAS? `ECE_TUNE_RVLAMSK*` and friends
   are indexed by vegetation type, not latitude, and every type leaks outside the
   boreal zone -- 38 % of deciduous needleleaf and 30 % of tundra by area lie
   equatorward of 50N. If those regions are unbiased or warm-biased, a lever
   justified by the boreal cold bias actively damages them. Types 4 and 9 are
   therefore split >50N vs <50N here.

Method. Period-clean: model `amip_presentday` (1990-2014) against ERA5 over the
same years, both on the model grid, land only. Each land cell is assigned the
single dominant surface type -- whichever of the high-vegetation type (weighted by
cvh), the low-vegetation type (cvl) or bare soil (1-cvl-cvh) has the largest
fraction -- so every cell is counted exactly once and the areas sum to the land
area. That is a coarser view than HTESSEL's own tiling, which lets several tiles
coexist in a cell, but it is the honest way to attribute a grid-box mean bias.

CAVEAT worth keeping in view: ERA5's land surface is HTESSEL, the same family as
ours. For 2m TEMPERATURE this matters much less than for the land-surface fields
(ERA5 assimilates screen-level observations, and T2m is strongly constrained by
them), but it is not a fully independent reference either.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF, ERA5_T2M

PD = list(range(1990, 2015))
SEAS = {'DJF': [11, 0, 1], 'MAM': [2, 3, 4], 'JJA': [5, 6, 7], 'SON': [8, 9, 10]}
NAMES = {0: 'bare soil', 1: 'Crops', 2: 'Short Grass', 3: 'Evgr Needleleaf',
         4: 'Dec Needleleaf', 5: 'Dec Broadleaf', 6: 'Evgr Broadleaf', 7: 'Tall Grass',
         8: 'Desert', 9: 'Tundra', 10: 'Irrig Crops', 11: 'Semidesert', 12: 'Ice Caps',
         13: 'Bogs/Marshes', 14: 'Inland Water', 15: 'Ocean', 16: 'Evgr Shrubs',
         17: 'Dec Shrubs', 18: 'Mixed Forest', 19: 'Interrupted Forest',
         20: 'Water/Land Mix'}

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def model_monthly(run, var, years):
    acc, n, lat, lon = None, 0, None, None
    for y in years:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        d = xr.open_dataset(f)
        a = d[var].values
        lat, lon = d[var].lat.values, d[var].lon.values
        d.close()
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, lat, lon) if n else (None, None, None)


def era5_on_grid(lat, lon):
    """ERA5 1990-2014 monthly climatology interpolated to the model grid."""
    d = xr.open_dataset(ERA5_T2M)
    a = d['t2m'].values
    la = d['latitude'].values if 'latitude' in d else d['lat'].values
    lo = d['longitude'].values if 'longitude' in d else d['lon'].values
    d.close()
    clim = a.reshape(len(a) // 12, 12, a.shape[1], a.shape[2]).mean(axis=0)
    da = xr.DataArray(clim, dims=('m', 'y', 'x'), coords={'y': la, 'x': lo})
    tlo = np.where(lon < 0, lon + 360, lon)
    return da.interp(y=xr.DataArray(np.clip(lat, la.min(), la.max()), dims='ny'),
                     x=xr.DataArray(tlo, dims='nx')).values


def dominant_type(lat, lon):
    """One surface type per land cell: the largest of cvh, cvl, bare."""
    def g(v):
        d = xr.open_dataset(f'{RT}/amip_pi_base/outdata/oifs/atm_remapped_1d_{v}_1d_1900-1900.nc')
        a = d[v].values[0]
        d.close()
        return a
    tvh, cvh, tvl, cvl = g('tvh'), g('cvh'), g('tvl'), g('cvl')
    bare = np.clip(1.0 - cvh - cvl, 0.0, 1.0)
    dom = np.where(cvh >= np.maximum(cvl, bare), np.round(tvh),
                   np.where(cvl >= bare, np.round(tvl), 0)).astype(int)
    return dom


t2m, lat, lon = model_monthly('amip_presentday', '2t', PD)
if t2m is None:
    raise SystemExit('amip_presentday 2t missing')
e5 = era5_on_grid(lat, lon)
dom = dominant_type(lat, lon)
w = np.cos(np.deg2rad(lat))[:, None] * np.ones_like(t2m[0])
land = lsm > 0.5
latg = np.broadcast_to(lat[:, None], t2m[0].shape)

bias = {s: (t2m[idx].mean(0) - e5[idx].mean(0)) for s, idx in SEAS.items()}


def row(mask, label, indent=''):
    if not mask.any():
        return
    a = w[mask].sum()
    vals = [np.average(bias[s][mask], weights=w[mask]) for s in SEAS]
    print(f'  {indent}{label:26s} {100*a/w[land].sum():6.2f} %' +
          ''.join(f'{v:>9.2f}' for v in vals))


print('Model amip_presentday MINUS ERA5, 1990-2014, land only, by dominant surface type.')
print('Negative = model too cold. Area is % of global land.\n')
print(f'  {"surface type":26s} {"area":>8s}' + ''.join(f'{s:>9s}' for s in SEAS))
print('  ' + '-' * 70)
order = sorted({int(t) for t in np.unique(dom[land])},
               key=lambda t: -w[land & (dom == t)].sum())
for t in order:
    m = land & (dom == t)
    if 100 * w[m].sum() / w[land].sum() < 0.4:
        continue
    row(m, f'{t:2d} {NAMES.get(t, "?")}')
print('  ' + '-' * 70)
row(land, 'ALL LAND')
row(land & (latg >= 50), 'all land >50N')
row(land & (latg < 50) & (latg > -50), 'all land 50S-50N')

print('\n\nThe leak test -- does a type\'s NON-BOREAL area share the bias?')
print('(these are the types a vegetation-indexed lever cannot separate)\n')
print(f'  {"":26s} {"area":>8s}' + ''.join(f'{s:>9s}' for s in SEAS))
print('  ' + '-' * 70)
for t in (3, 4, 9, 16, 17, 13):
    m = land & (dom == t)
    if not m.any():
        continue
    print(f'  {t:2d} {NAMES[t]}')
    row(m & (latg >= 50), 'poleward of 50N', '   ')
    row(m & (latg < 50), 'equatorward of 50N', '   ')
print('\n  A type whose equatorward area is UNBIASED or WARM cannot safely be given a')
print('  lever justified by the boreal cold bias -- the same namelist entry acts on both.')

# ---------------------------------------------------------------------------
# RESULT (2026-08-04). Two findings, the second much larger than the question
# this script was written to answer.
#
# 1. THE LEAK IS BENIGN. The non-boreal areas of types 4 and 9 are ALSO cold
#    biased (equatorward of 50N: type 4 JJA -0.68, type 9 JJA -1.06), so a
#    warming lever applied to those types is directionally correct everywhere it
#    acts. The 30-38 % leak does not damage the regions it reaches.
#
# 2. THE COLD BIAS IS NOT BOREAL. Every one of the 20 surface types is cold in
#    every season -- not one is warm. Global land is DJF -1.73, MAM -1.60,
#    JJA -1.35, SON -1.39. The Siberian JJA bias (-2.58) is about TWICE the
#    global-land JJA bias, so there is a real boreal enhancement, but most of the
#    Siberian bias is shared with all land.
#
#    Rough decomposition of the -2.6 K Siberian JJA bias:
#        ~0.7 K  present even over PRESCRIBED-SST ocean (-0.72 K there)
#        ~0.7 K  additional over land generally
#        ~1.2 K  additional over the boreal zone
#    G4 has recovered ~0.95 K, i.e. most of the boreal-specific budget -- which
#    is very likely WHY the F-series saturated (F5 = F4 to within noise). A
#    boreal-indexed lever cannot reach the other ~1.4 K.
#
# CAVEATS. ERA5's land surface is HTESSEL, though 2m temperature is strongly
# constrained by assimilated screen observations, so it is a fairer reference
# here than for the land-surface fields. Orographic differences between the
# 0.25-degree ERA5 and the model grid contribute some of the land bias. Over
# ocean, part of the -0.72 K may be an SST-dataset difference rather than an air
# temperature bias -- the model's own air-sea difference is -1.53 K globally,
# at the cold edge of the observed -0.5 to -1.5 K range but not outside it.
# ---------------------------------------------------------------------------
