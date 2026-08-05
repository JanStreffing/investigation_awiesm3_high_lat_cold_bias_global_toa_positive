"""Is the -22 K winter soil bias of the I-series real, and was there a bias before it?

WHY THIS EXISTS. `monthly` diagnostics showed the round-15 snow-depletion scheme drives
Siberian soil layer 1 to -23.8 K against the control by January, with the divergence starting
in October when snow is shallow and the tanh cover formula bites hardest. That was reported
as a DELTA, which says nothing about whether it is a NEW error or the removal of an old one.
This campaign has already made exactly that mistake once -- the retracted claim that excess
snow cover was "propping DJF up ~2.7 K" -- so the control must be scored against a reference
before the delta means anything.

TWO QUESTIONS, ANSWERED SEPARATELY:
  1. Is the CONTROL's winter soil temperature right?  (was there a bias before?)
  2. Does the I-series move toward or away from the reference?

HOW WELL CONSTRAINED IS ERA5 HERE?  This is the crux, and the answer is: much less than for
T2m, so ERA5 alone cannot settle it.
  * ERA5's land surface IS HTESSEL -- the same model family, the same snow scheme, the same
    soil thermal conductivity formulation. It is a sibling, not an independent witness.
  * ERA5 assimilates SCREEN-LEVEL observations (2 m T and RH) into the soil through a
    simplified Extended Kalman Filter, and that analysis targets soil MOISTURE far more
    strongly than soil temperature. Deep soil temperature is essentially model-determined.
  * There is no dense assimilated soil-temperature observing network in Siberia.
So agreement with ERA5 is weak evidence, but DISAGREEMENT is still informative: if our
control already differs from a sibling HTESSEL that is nudged by screen observations, that is
a real problem, and if the I-series moves further away it is very unlikely to be a correction.

THE INDEPENDENT PHYSICAL CHECK, which does not depend on ERA5 at all. The winter air-soil
temperature difference under snow -- the "thermal offset" or nival offset -- is a measured
quantity in the permafrost literature, from Russian station networks (RIHMI-WDC, Bulygina
et al.) and GTN-P boreholes. Under a deep Siberian snowpack it is characteristically
+10 to +20 K (soil warmer than air), and it is NEVER ~0: a metre-plus of snow is an excellent
insulator and the ground simply does not track air temperature. So:
      offset ~ 0 K  => the model has effectively no snow insulation, which is unphysical
      offset >> 20 K => over-insulated
That bound is what makes this diagnostic decisive regardless of ERA5's limitations.

PERIOD-CLEAN: model `amip_presentday` (1990-2014) against ERA5 over the same years. The
tuning runs are 1872-1915 and are compared to each other, never directly to ERA5.
"""
import numpy as np, xarray as xr, os, subprocess, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF

PD = list(range(1990, 2015))
DJF = [11, 0, 1]
WORK = '/work/ab0246/a270092/obs/era5/soil'
CODES = {'stl1': 139, 'stl2': 170, 'stl3': 183, 'stl4': 236}
TUNING = [('control', 'amip_pi_base'), ('G4', 'amip_G4_tundra'),
          ('I1', 'amip_I1_scf'), ('I3', 'amip_I3_scf_sdor')]

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def model_clim(run, var, years):
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


def era5_clim(var, lat, lon):
    """ERA5 monthly climatology for `var` over PD, on the model grid."""
    os.makedirs(WORK, exist_ok=True)
    out = f'{WORK}/e5_{var}_clim.nc'
    if not os.path.exists(out):
        src = [f'/pool/data/ERA5/E5/sf/an/1M/{CODES[var]}/E5sf00_1M_{y}_{CODES[var]}.grb'
               for y in PD]
        src = [s for s in src if os.path.exists(s)]
        if not src:
            return None
        cmd = ('module load cdo >/dev/null 2>&1; '
               f'cdo -s -f nc -ymonmean -setgridtype,regular -cat "{" ".join(src)}" {out}')
        subprocess.run(['bash', '-lc', cmd], check=False)
    if not os.path.exists(out):
        return None
    d = xr.open_dataset(out)
    # cdo writes (time, depth, lat, lon) for soil fields -- squeeze the singleton depth
    v = [k for k in d.data_vars if d[k].ndim >= 3 and 'bnds' not in k][0]
    a = d[v].squeeze()
    la = d['lat'].values if 'lat' in d else d['latitude'].values
    lo = d['lon'].values if 'lon' in d else d['longitude'].values
    da = xr.DataArray(a.values, dims=('m', 'y', 'x'), coords={'y': la, 'x': lo})
    if la[0] > la[-1]:
        da = da.isel(y=slice(None, None, -1))
    tlo = np.where(lon < 0, lon + 360, lon)
    r = da.interp(y=xr.DataArray(np.clip(lat, min(la), max(la)), dims='ny'),
                  x=xr.DataArray(tlo, dims='nx')).values
    d.close()
    return r


t2m, lat, lon = model_clim('amip_presentday', '2t', PD)
if t2m is None:
    raise SystemExit('amip_presentday 2t missing')
la2 = np.broadcast_to(lat[:, None], t2m[0].shape)
box = (la2 >= 55) & (la2 <= 75) & \
      np.broadcast_to(((lon >= 60) & (lon <= 180))[None, :], t2m[0].shape) & (lsm > 0.5)
w = np.cos(np.deg2rad(lat))[:, None] * np.ones_like(t2m[0])


def bm(a, m=box):
    return np.average(a[m], weights=w[m])


print(__doc__.split('PERIOD-CLEAN')[0])
print('=' * 78)
print('QUESTION 1 -- was the control already biased?  Siberian land, DJF, 1990-2014.')
print('=' * 78)
e5_t2 = era5_clim('stl1', lat, lon)  # warm the cache path; t2m handled below
rows = []
for v in ('stl1', 'stl2', 'stl3', 'stl4'):
    mv, _, _ = model_clim('amip_presentday', v, PD)
    if mv is None:
        print(f'  {v}: model output missing')
        continue
    ev = era5_clim(v, lat, lon)
    m_ = bm(mv[DJF].mean(0))
    if ev is None:
        print(f'  {v}: model {m_-273.15:7.2f} C   ERA5 unavailable')
        continue
    e_ = bm(ev[DJF].mean(0))
    rows.append((v, m_, e_))
    print(f'  {v}: model {m_-273.15:7.2f} C   ERA5 {e_-273.15:7.2f} C   bias {m_-e_:+6.2f} K')

m_t2 = bm(t2m[DJF].mean(0))
print(f'\n  2m air (model, same box/period): {m_t2-273.15:7.2f} C')
if rows:
    print(f'  => control THERMAL OFFSET (soil L1 - air) = {rows[0][1]-m_t2:+.2f} K')
    print(f'     ERA5 thermal offset                    = {rows[0][2]-m_t2:+.2f} K'
          '   [same air, so this isolates the soil]')

print()
print('=' * 78)
print('QUESTION 2 -- does the I-series move toward or away?  1872-1915, vs control.')
print('=' * 78)
res = {}
for nm, r in TUNING:
    mv, _, _ = model_clim(r, 'stl1', range(1872, 1916))
    ta, _, _ = model_clim(r, '2t', range(1872, 1916))
    if mv is None:
        continue
    res[nm] = (bm(mv[DJF].mean(0)), bm(ta[DJF].mean(0)))
print(f'  {"run":10s}{"soil L1 [C]":>13s}{"2m air [C]":>12s}{"offset [K]":>12s}')
for nm, _ in TUNING:
    if nm not in res:
        continue
    s, a = res[nm]
    print(f'  {nm:10s}{s-273.15:13.2f}{a-273.15:12.2f}{s-a:12.2f}')
print("""
  READING IT. The thermal offset is the physically bounded quantity: Russian station and
  GTN-P borehole data put the winter soil-minus-air difference under a deep Siberian
  snowpack at roughly +10 to +20 K, and never near zero. An offset that collapses toward
  0 K means the model has lost its snow insulation altogether, which no observation
  supports -- so a large NEGATIVE soil change is a new error, not the removal of an old one.
  Judge the I-series on the OFFSET, not on the soil temperature alone, because the offset
  is insensitive to whatever air-temperature bias the run also carries.
""")
