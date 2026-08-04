"""The cold bias through the whole column: soil layer 4 -> skin -> 2 m -> 100 hPa.

Why this exists. `bias_by_tile.py` found the model is cold over EVERY land surface
type in every season, and that ~0.7 K of it survives over near-sea-level land AND
over prescribed-SST ocean -- i.e. a component that is neither boreal, nor
elevation-driven, nor a land-surface parameter. The obvious next question is where
in the vertical that sits, because the answer discriminates sharply between
candidate causes:

  * bias confined to soil + skin + 2 m, vanishing by 850 hPa
        -> a SURFACE or DIAGNOSTIC problem (the 2 m interpolation, the skin layer,
           soil heat storage). Fixable at the surface.
  * bias roughly uniform through the troposphere
        -> an ATMOSPHERIC problem (radiation, lapse rate, advection). No surface
           lever will touch it, and the whole campaign's premise would need
           revisiting.
  * bias growing downward from a small free-tropospheric value
        -> boundary-layer coupling: the surface is decoupled from an atmosphere
           that is roughly right.

The column is assembled from, in order:
    stl4, stl3, stl2, stl1   soil temperature, layers 4->1 (2.89 m -> 7 cm)
    skt                      skin temperature
    2t                       2 m temperature
    t on pressure levels     1000, 925, 850, 700, 500, 300, 200, 100 hPa

Period-clean throughout: model `amip_presentday` (1990-2014) vs ERA5 over the same
years. Reported for global land, and separately for LOW-LYING land (<500 m) so the
orographic-mismatch artefact identified in `bias_by_tile.py` is bounded, and for
the Siberian box for continuity with the rest of the campaign.

CAVEAT. ERA5's soil and skin fields are HTESSEL output, the same scheme family as
ours, so the soil rows are model-vs-sibling and are NOT independent evidence --
they show whether OUR soil differs from a well-initialised version of the same
scheme, which is still useful but is not observation. 2 m temperature and the
pressure-level temperatures are strongly constrained by assimilated observations
and are the trustworthy rows.
"""
import numpy as np, xarray as xr, os, subprocess, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

PD = list(range(1990, 2015))
SEAS = {'DJF': [11, 0, 1], 'JJA': [5, 6, 7]}
W = '/work/ab0246/a270092/obs/era5/column'
POOL_SF = '/pool/data/ERA5/E5/sf/an/1M'
POOL_PL = '/pool/data/ERA5/E5/pl/an/1M'
# ERA5 surface params: soil layers 1-4, skin, 2 m
SF = [('139', 'stl1'), ('170', 'stl2'), ('183', 'stl3'), ('236', 'stl4'),
      ('235', 'skt'), ('167', '2t')]
PLEV = [1000, 925, 850, 700, 500, 300, 200, 100]

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def model_clim(var, pl=False):
    acc, n, lat, lon = None, 0, None, None
    for y in PD:
        tag = 'pl_t_1m_pl' if pl else f'{var}_1m'
        f = f'{RT}/amip_presentday/outdata/oifs/atm_remapped_1m_{tag}_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        d = xr.open_dataset(f)
        key = 't' if pl else var
        a = d[key].values
        lat, lon = d[key].lat.values, d[key].lon.values
        lev = (d[key].pressure_levels.values if pl and 'pressure_levels' in d[key].coords
               else (d[key].plev.values if pl and 'plev' in d[key].coords else None))
        d.close()
        acc = a if acc is None else acc + a
        n += 1
    if not n:
        return None, None, None, None
    return acc / n, lat, lon, (lev if pl else None)


def build_era5():
    """One-off: ERA5 climatologies on the model grid."""
    os.makedirs(W, exist_ok=True)
    grid = os.path.join(W, 'modelgrid.txt')
    if not os.path.exists(grid):
        src = f'{RT}/amip_presentday/outdata/oifs/atm_remapped_1m_2t_1m_1990-1990.nc'
        subprocess.run(f'cdo -s griddes {src} > {grid}', shell=True, check=True)
    for code, name in SF:
        out = os.path.join(W, f'e5_{name}.nc')
        if os.path.exists(out):
            continue
        files = ' '.join(f'{POOL_SF}/{code}/E5sf00_1M_{y}_{code}.grb' for y in PD
                         if os.path.exists(f'{POOL_SF}/{code}/E5sf00_1M_{y}_{code}.grb'))
        if not files:
            print(f'  !! ERA5 {code} ({name}) has no input'); continue
        subprocess.run(f'cdo -s -O -f nc -remapbil,{grid} -ymonmean -setgridtype,regular '
                       f'-cat "{files}" {out}', shell=True)
    out = os.path.join(W, 'e5_t_pl.nc')
    if not os.path.exists(out):
        lev = ','.join(str(p * 100) for p in PLEV)
        files = ' '.join(f'{POOL_PL}/130/E5pl00_1M_{y}_130.grb' for y in PD
                         if os.path.exists(f'{POOL_PL}/130/E5pl00_1M_{y}_130.grb'))
        if files:
            subprocess.run(f'cdo -s -O -f nc -remapbil,{grid} -ymonmean -sellevel,{lev} '
                           f'-setgridtype,regular -cat "{files}" {out}', shell=True)


def e5(name):
    p = os.path.join(W, f'e5_{name}.nc')
    if not os.path.exists(p):
        return None
    d = xr.open_dataset(p)
    # pick the real field, not a time/bounds variable (cdo emits several)
    cand = [k for k in d.data_vars
            if d[k].ndim >= 3 and not np.issubdtype(d[k].dtype, np.datetime64)]
    if not cand:
        d.close(); return None
    a = np.squeeze(d[cand[0]].values)
    la = d['lat'].values
    d.close()
    if la[0] > la[-1]:
        a = a[:, ::-1] if a.ndim == 3 else a[:, :, ::-1]
    return a


print('Building ERA5 column climatologies (first run only, a few minutes)...')
build_era5()

m2t, lat, lon, _ = model_clim('2t')
if m2t is None:
    raise SystemExit('model amip_presentday missing')
land = lsm > 0.5
w = np.cos(np.deg2rad(lat))[:, None] * np.ones_like(m2t[0])

# orography for the low-lying subset
oro = None
op = '/work/ab0246/a270092/obs/era5/column/orog_rg.nc'
if os.path.exists(op):
    d = xr.open_dataset(op)
    _c = [k for k in d.data_vars if d[k].ndim >= 2 and not np.issubdtype(d[k].dtype, np.datetime64)]
    oro = np.squeeze(d[_c[0]].values) / 9.80665; d.close()
    if oro.ndim == 3: oro = oro[0]

BOX = ((55, 75), (60, 180))
ys = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
l180 = ((lon + 180) % 360) - 180
xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
sib = np.zeros_like(land); sib[np.ix_(np.where(ys)[0], np.where(xs)[0])] = True
sib &= land

masks = [('global land', land)]
if oro is not None:
    masks.append(('land <500 m', land & (oro < 500)))
masks.append(('Siberian box', sib))

rows = []
for name in ('stl4', 'stl3', 'stl2', 'stl1', 'skt', '2t'):
    mm, _, _, _ = model_clim(name)
    ee = e5(name)
    if mm is None or ee is None:
        rows.append((name, None)); continue
    rows.append((name, (mm, ee)))

mpl, _, _, lev = model_clim('t', pl=True)
epl = e5('t_pl')
# ERA5 was written with its own level ordering; align it to the model's explicitly
elev = None
_p = os.path.join(W, 'e5_t_pl.nc')
if os.path.exists(_p):
    _d = xr.open_dataset(_p); elev = _d['plev'].values; _d.close()

print('\nMODEL amip_presentday MINUS ERA5, 1990-2014 [K].  Negative = model too cold.\n')
for sname, sidx in SEAS.items():
    print(f'  === {sname} ===')
    print(f'    {"level":12s}' + ''.join(f'{n:>16s}' for n, _ in masks))
    if mpl is not None and epl is not None and lev is not None:
        order = np.argsort(-np.asarray(lev))
        for k in order:
            p = lev[k] / 100.0
            if int(round(p)) not in PLEV or elev is None:
                continue
            ke = int(np.argmin(np.abs(np.asarray(elev) - lev[k])))   # match by VALUE, not index
            if abs(elev[ke] - lev[k]) > 1.0:
                continue
            b = mpl[sidx, k].mean(0) - epl[sidx, ke].mean(0)
            print(f'    {int(round(p)):5d} hPa   ' +
                  ''.join(f'{np.average(b[m], weights=w[m]):16.2f}' for _, m in masks))
    for name, dat in reversed(rows):
        if dat is None:
            print(f'    {name:12s}' + ''.join(f'{"--":>16s}' for _ in masks)); continue
        mm, ee = dat
        b = mm[sidx].mean(0) - ee[sidx].mean(0)
        lbl = {'2t': '2 m', 'skt': 'skin', 'stl1': 'soil L1 7cm',
               'stl2': 'soil L2 28cm', 'stl3': 'soil L3 1m', 'stl4': 'soil L4 2.9m'}[name]
        print(f'    {lbl:12s}' + ''.join(f'{np.average(b[m], weights=w[m]):16.2f}' for _, m in masks))
    print()
print('  Reading it: a bias confined to soil/skin/2m is a SURFACE or DIAGNOSTIC problem;')
print('  one that persists to 500 hPa is ATMOSPHERIC and no surface lever will fix it.')
print('  ERA5 soil/skin rows are HTESSEL-vs-HTESSEL and are not independent evidence.')

# ---------------------------------------------------------------------------
# RESULT (2026-08-04). This overturns the campaign's working premise.
#
# THE BIAS IS NOT AT THE SURFACE. Reading the JJA column over low-lying land
# (where the orographic artefact is bounded):
#     soil L1 -0.36, L2 -0.38, L3 -0.72, L4 -0.89
#     skin    -0.53
#     2 m     -1.00
#     850 hPa -0.92,  700 -0.89,  500 -1.24,  300 -1.85,  200 -2.97,  100 -1.70
#
# The soil and skin are the LEAST biased parts of the column. The free
# troposphere is cold by ~0.9 K, and 2 m sits only ~0.1-0.4 K below it. So the
# near-surface cold is largely INHERITED from a cold atmosphere rather than
# generated at the surface -- and no surface parameter can fix that.
#
# The upper troposphere is much worse: -2.8 to -3.0 K at 200 hPa globally and
# -4.9 K over Siberia in JJA. A tropopause-region cold bias of that size is a
# known IFS/EC-Earth trait and is a separate problem from the surface work.
#
# WHAT REMAINS SURFACE-CONFINED. Siberia JJA: 2 m -2.04, 925 hPa -2.17,
# 850 -1.99, but 700 hPa only -1.22. So the Siberian boundary layer carries
# ~0.8 K on top of the ~1.2 K tropospheric bias. G4 recovered ~0.95 K, i.e.
# approximately the whole boundary-layer-confined budget -- which is a much better
# explanation of the F-series saturation (F5 = F4) than any parameter limit.
#
# CONSEQUENCE. "Warm the land surface by 0.7 K" is not a surface problem. The
# atmosphere above the land is cold by ~0.9 K in the lower free troposphere, and
# the surface is in near-equilibrium with it. Surface levers can only reach the
# ~0.2-0.4 K by which 2 m departs from 850 hPa.
#
# CAVEATS. The 1000 hPa row is unreliable over land -- much of that level is
# below ground and both datasets extrapolate; use 925/850 as the lowest honest
# levels. ERA5 soil and skin are HTESSEL output, so those rows are
# model-vs-sibling; the pressure-level temperatures and 2 m are observation-
# constrained and are the rows to trust. Note also the sign puzzle worth keeping
# in view: global net TOA is +0.52 W/m2 (energy coming IN) while the troposphere
# is cold -- in AMIP the prescribed SST absorbs that imbalance without warming, so
# a positive TOA imbalance and a cold atmosphere are not contradictory here.
# ---------------------------------------------------------------------------
