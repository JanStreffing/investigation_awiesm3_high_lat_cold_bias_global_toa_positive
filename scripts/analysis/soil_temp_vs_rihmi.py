"""The winter soil thermal offset against DIRECT Russian station observations.

THE QUESTION. The round-15 snow-depletion scheme (I1, I3) drives Siberian soil layer 1 to
-30.5 C in DJF against a control of -8.2 C. Reported as a delta that is meaningless: it says
nothing about whether it is a NEW error or the removal of an old one. `soil_temp_vs_obs.py`
answered that against ERA5 -- control biased by only +1.08 K, I-series thermal offset
collapsing from +21 K to +1.4 K -- but ERA5's land surface IS HTESSEL, the same model family
with the same snow scheme, so it is a sibling rather than a witness.

THIS SCRIPT USES REAL OBSERVATIONS. RIHMI-WDC (Sherstyukov v3) daily soil temperature at
Russian meteorological stations, 1963-2024, at 12 depths from 2 to 320 cm, measured UNDER
NATURAL COVER -- the snowpack is left intact, which is exactly the condition whose insulation
is in question. 43 of the 110 stations fall in the campaign's Siberian box (55-75N, 60-180E).

WHAT THE ARCHIVE CAN AND CANNOT SETTLE.
  * Depths 2/5/10/15 cm carry NO winter data: the Russian network withdraws the shallow
    Savinov thermometers for the cold season. Winter observations begin at 20 cm. So model
    `stl1` (0-7 cm), which is where the I-series damage was reported, has NO direct
    observational counterpart. `stl2` (7-28 cm) is validated instead, and it is tightly
    coupled to `stl1`, so it bounds it.
  * Observation time is not uniform (80-320 cm read once daily near 14 h LST; 20 and 40 cm
    follow synoptic hours in the warm season). For a DJF mean at 20 cm and below this is a
    minor error compared with the effect being tested.
  * Quality flags are used strictly: only `tsoil_qc == 0` ("value is reliable"). The archive
    authors deliberately leave suspect values in place for the user to screen.

METHOD -- CO-LOCATED, not box-averaged. Model and ERA5 are sampled AT THE STATION POINTS by
nearest-neighbour, so the comparison is not contaminated by the box's forest/bog/tundra mix
being different from where the stations happen to sit. That is the main representativeness
control available without a land-cover weighting.

THE DIAGNOSTIC IS THE THERMAL OFFSET, soil minus 2 m air, because:
  (a) it is insensitive to whatever air-temperature bias a run carries -- I1's air is itself
      2.8 K colder than the control's, which would confound a raw soil comparison;
  (b) it is the direct measure of snow insulation, which is the physics under test;
  (c) each dataset uses its OWN air temperature, so no cross-dataset assumption enters.
Observed offset = station soil minus ERA5 T2m (ERA5 assimilates screen observations densely,
so T2m is the one field where it IS well constrained). Model offset = model soil minus model
T2m at the same points.

PERIOD-CLEAN: `amip_presentday` (1990-2014) against stations over the same years. The tuning
runs are 1872-1915 and are compared on OFFSET only, never on absolute temperature.
"""
import numpy as np, xarray as xr, pandas as pd, os, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF, ERA5_T2M

OBSF = '/work/ab0246/a270092/obs/RIHMI-WDC/RIHMI-WDC_soil_temperature_v3_1963-2024.nc'
PD = list(range(1990, 2015))
DJF = [11, 0, 1]
BOX = (55, 75, 60, 180)
# model layer -> (depth range, observation depths that bracket its centre)
LAYERS = {'stl2': ((0.07, 0.28), [0.20]),
          'stl3': ((0.28, 1.00), [0.40, 0.80]),
          'stl4': ((1.00, 2.89), [1.20, 1.60, 2.40])}
TUNING = [('control', 'amip_pi_base'), ('G4', 'amip_G4_tundra'),
          ('I1', 'amip_I1_scf'), ('I3', 'amip_I3_scf_sdor')]


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


def at_stations(field, lat, lon, slat, slon):
    """Nearest-neighbour sample of a (12, ny, nx) field at station points."""
    lo = np.where(lon < 0, lon + 360, lon)
    order = np.argsort(lo)
    out = np.empty((field.shape[0], len(slat)))
    for k, (la, ln) in enumerate(zip(slat, slon)):
        j = order[np.clip(np.searchsorted(lo[order], ln % 360), 0, len(lo) - 1)]
        i = int(np.argmin(np.abs(lat - la)))
        out[:, k] = field[:, i, j]
    return out


# ---- observations -----------------------------------------------------------
d = xr.open_dataset(OBSF)
slat, slon = d.lat.values, d.lon.values
inbox = (slat >= BOX[0]) & (slat <= BOX[1]) & (slon >= BOX[2]) & (slon <= BOX[3])
tt = pd.DatetimeIndex(d.time.values)
sel = (tt.year >= 1990) & (tt.year <= 2014) & tt.month.isin([12, 1, 2])
ts = np.where(d.tsoil_qc.values == 0, d.tsoil.values, np.nan)   # strict QC
dep = d.depth.values
obs = {}
for k, dd in enumerate(dep):
    v = ts[inbox][:, sel, k]
    if np.isfinite(v).sum() > 1000:
        obs[round(float(dd), 2)] = np.nanmean(v, axis=1)         # per-station DJF mean
slat_b, slon_b = slat[inbox], slon[inbox]
d.close()

# ---- ERA5 T2m at the same stations ------------------------------------------
e = xr.open_dataset(ERA5_T2M)
ea = e['t2m'].values
ela = e['latitude'].values if 'latitude' in e else e['lat'].values
elo = e['longitude'].values if 'longitude' in e else e['lon'].values
e.close()
eclim = ea.reshape(len(ea) // 12, 12, ea.shape[1], ea.shape[2]).mean(axis=0)
if ela[0] > ela[-1]:
    eclim = eclim[:, ::-1, :]
    ela = ela[::-1]
e_t2 = at_stations(eclim, ela, elo, slat_b, slon_b)[DJF].mean(0) - 273.15

print(__doc__.split('PERIOD-CLEAN')[0])
print('=' * 80)
print(f'Siberian box {BOX}: {inbox.sum()} RIHMI stations, DJF 1990-2014, qc==0 only.')
print('=' * 80)

# ---- model at the same stations ---------------------------------------------
t2m, mlat, mlon = model_clim('amip_presentday', '2t', PD)
if t2m is None:
    raise SystemExit('amip_presentday 2t missing')
m_t2 = at_stations(t2m, mlat, mlon, slat_b, slon_b)[DJF].mean(0) - 273.15

print(f'\n  2 m air at the stations:  ERA5 {np.nanmean(e_t2):7.2f} C'
      f'   model {np.nanmean(m_t2):7.2f} C   bias {np.nanmean(m_t2-e_t2):+.2f} K')
print('\n  Soil temperature and THERMAL OFFSET (soil - own 2 m air):\n')
print(f'  {"layer":6s} {"obs depth":>10s} {"OBS T":>8s} {"model T":>9s} {"bias":>7s}'
      f' | {"OBS offset":>11s} {"model offset":>13s} {"error":>8s}')
print('  ' + '-' * 88)
for lay, (rng, depths) in LAYERS.items():
    mv, _, _ = model_clim('amip_presentday', lay, PD)
    if mv is None:
        continue
    m_s = at_stations(mv, mlat, mlon, slat_b, slon_b)[DJF].mean(0) - 273.15
    have = [dd for dd in depths if dd in obs]
    if not have:
        continue
    o_s = np.nanmean(np.vstack([obs[dd] for dd in have]), axis=0)
    ok = np.isfinite(o_s) & np.isfinite(m_s)
    o_off = np.nanmean(o_s[ok] - e_t2[ok])
    m_off = np.nanmean(m_s[ok] - m_t2[ok])
    print(f'  {lay:6s} {"/".join(f"{x*100:.0f}" for x in have)+" cm":>10s}'
          f' {np.nanmean(o_s[ok]):8.2f} {np.nanmean(m_s[ok]):9.2f}'
          f' {np.nanmean(m_s[ok]-o_s[ok]):+7.2f} | {o_off:+11.2f} {m_off:+13.2f}'
          f' {m_off-o_off:+8.2f}')

print(f'\n  stl1 (0-7 cm): NO winter observations exist -- shallow thermometers are withdrawn')
print(f'  for the cold season. stl2 above is the shallowest verifiable layer.')

# ---- the tuning runs, on OFFSET only ----------------------------------------
print()
print('=' * 80)
print('THE TUNING RUNS (1872-1915). Offset only -- absolutes are not period-comparable.')
print('=' * 80)
print(f'\n  {"run":9s}{"stl1 off":>10s}{"stl2 off":>10s}{"stl3 off":>10s}{"stl4 off":>10s}')
for nm, r in TUNING:
    ta, _, _ = model_clim(r, '2t', range(1872, 1916))
    if ta is None:
        continue
    a = at_stations(ta, mlat, mlon, slat_b, slon_b)[DJF].mean(0)
    row = f'  {nm:9s}'
    for lay in ('stl1', 'stl2', 'stl3', 'stl4'):
        sv, _, _ = model_clim(r, lay, range(1872, 1916))
        if sv is None:
            row += f'{"--":>10s}'
            continue
        s = at_stations(sv, mlat, mlon, slat_b, slon_b)[DJF].mean(0)
        row += f'{np.nanmean(s-a):10.2f}'
    print(row)
print(f'\n  OBSERVED offsets for comparison:', end='')
for lay, (rng, depths) in LAYERS.items():
    have = [dd for dd in depths if dd in obs]
    if not have:
        continue
    o_s = np.nanmean(np.vstack([obs[dd] for dd in have]), axis=0)
    ok = np.isfinite(o_s)
    print(f'  {lay} {np.nanmean(o_s[ok]-e_t2[ok]):+.1f}', end='')
print("""

  READING IT. The offset is snow insulation made visible, and it is now measured, not
  assumed. A model offset near ZERO means the soil tracks the air, i.e. the snowpack has
  stopped insulating -- which the station record directly refutes. Judge the I-series
  against the OBSERVED offset column, not against the control.
""")
