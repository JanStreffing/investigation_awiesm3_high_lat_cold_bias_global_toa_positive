"""Where does the ocean heat uptake enter, and where is it stored?

THE HYPOTHESIS BEING TESTED, stated before looking.  The arms bury ~0.8 W/m2 in the ocean
while T2m stays flat.  A sustained uptake cannot sit in a tropical or subtropical mixed layer,
because that layer re-equilibrates with the atmosphere in a few years and simply radiates the
excess back.  Sustained storage requires the heat to be REMOVED from contact with the surface,
which happens where deep and intermediate water forms: the Southern Ocean, the subpolar North
Atlantic and the Nordic Seas.  Those are the three boxes this campaign has been tuning
against by surface shortwave RMSE.  If the hypothesis holds, the storage should concentrate
at high latitude and at depth, and the tropics should show little NET accumulation despite
being where the surface flux is largest.

TWO DIFFERENT QUESTIONS, kept separate here because conflating them is the easy mistake.
  ENTRY   is the net surface heat flux, FESOM `fh`, by latitude band.  It says where energy
          crosses the air-sea interface, but circulation moves it afterwards, so a band with
          large entry need not store anything.
  STORAGE is the change in heat content, from 3D `temp`, by band AND by depth.  This is the
          one that answers the question, because it is where the energy actually ends up.

METHOD.  Heat content is rho*cp*sum(T * nod_area[k,i] * dz[k]) over nodes i and levels k,
with nod_area taken per level so shoaling bathymetry is handled, and dz from the zbar
interfaces.  Storage is the late window minus the early window, divided by the elapsed time
and the WHOLE planet's area, so every number is directly comparable to the TOA imbalance in
W/m2.  Windows are 8 years each to average out interannual noise.

TRAP.  nod_area has nz levels of interfaces; layer k lies between interface k and k+1.

SIGN OF `fh`.  The file carries no sign convention beyond long_name "surface heat flux", and
as stored it gives the tropics -13.7 W/m2 and the Southern Ocean +18.6.  The tropical ocean
gains heat at the surface in the annual mean and the Southern Ocean loses it, so as stored
`fh` is positive UPWARD, out of the ocean.  It is negated below.  With that sign the ocean
gains +1.16 W/m2 of global area, which sits alongside the storage and TOA figures instead of
contradicting them; taken as written it would have implied the ocean was losing heat while
its heat content rose.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESH = '/work/ab0246/a270092/input/fesom2/core3'
RHO, CP = 1027.0, 3996.0
A_EARTH = 5.101e14
SEC_YR = 365.25 * 86400.0

BANDS = [('90-60S  Antarctic/SO', -90, -60),
         ('60-45S  subantarctic', -60, -45),
         ('45-30S', -45, -30),
         ('30S-30N  tropics', -30, 30),
         ('30-45N', 30, 45),
         ('45-60N  subpolar NA', 45, 60),
         ('60-90N  Arctic/Nordic', 60, 90)]
DEPTHS = [('0-100 m', 0, 100), ('100-700 m', 100, 700),
          ('700-2000 m', 700, 2000), ('>2000 m', 2000, 1e9)]

ARM = os.environ.get('ARM', '11Q')
EARLY = list(range(1352, 1360))
LATE = list(range(1382, 1390))

root = f'{R092}/{ARM}' if os.path.isdir(f'{R092}/{ARM}') else \
       [p for p in glob.glob(f'{R092}/*{ARM}*') if os.path.isdir(p)][0]

with xr.open_dataset(f'{MESH}/fesom.mesh.diag.nc', decode_times=False) as m:
    lat = m['lat'].values
    if np.abs(lat).max() < 4:
        lat = np.rad2deg(lat)
    nod_area = m['nod_area'].values          # (nz, nod2)
    zbar = m['nz'].values                    # (nz,) interfaces, negative downward
z = np.abs(zbar)
dz = np.diff(z)                              # (nz-1,) layer thicknesses
nlev = len(dz)
area = nod_area[:nlev, :]                    # (nlev, nod2)
vol = area * dz[:, None]                     # (nlev, nod2) m3
zmid = 0.5 * (z[:-1] + z[1:])
print(f'arm {ARM}   nodes {len(lat)}   layers {nlev}   '
      f'ocean volume {vol.sum():.4e} m3   depth to {z[-1]:.0f} m')


def annual_temp(y):
    f = glob.glob(f'{root}/outdata/fesom/temp.fesom.{y}.nc')
    if not f:
        return None
    with xr.open_dataset(f[0], decode_times=False) as d:
        n = 'temp' if 'temp' in d.data_vars else list(d.data_vars)[-1]
        return d[n].mean(dim='time').values.T.astype(np.float64)   # -> (nz, nod2)


def window_mean(years):
    acc, n = None, 0
    for y in years:
        a = annual_temp(y)
        if a is None:
            continue
        a = a[:nlev, :]
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, n) if n else (None, 0)


T0, n0 = window_mean(EARLY)
T1, n1 = window_mean(LATE)
if T0 is None or T1 is None:
    raise SystemExit(f'  insufficient temp output ({n0} early, {n1} late)')
dt_yr = np.mean(LATE) - np.mean(EARLY)
print(f'  windows {EARLY[0]}-{EARLY[-1]} ({n0} yr) -> {LATE[0]}-{LATE[-1]} ({n1} yr), '
      f'{dt_yr:.0f} yr apart\n')

dT = np.where(np.isfinite(T1) & np.isfinite(T0), T1 - T0, 0.0)
dH = RHO * CP * dT * vol                      # J, per (level, node)
to_wm2 = 1.0 / (dt_yr * SEC_YR * A_EARTH)

print('  STORAGE: change in ocean heat content, W/m2 of GLOBAL area\n')
hdr = f'  {"band":24s}' + ''.join(f'{d[0]:>12s}' for d in DEPTHS) + f'{"total":>10s}'
print(hdr); print('  ' + '-' * (len(hdr) - 2))
band_tot = []
for name, la0, la1 in BANDS:
    sel = (lat >= la0) & (lat < la1)
    row = []
    for _dn, d0, d1 in DEPTHS:
        lk = (zmid >= d0) & (zmid < d1)
        row.append(dH[np.ix_(lk, sel)].sum() * to_wm2)
    band_tot.append((name, row, sum(row)))
    print(f'  {name:24s}' + ''.join(f'{v:+12.4f}' for v in row) + f'{sum(row):+10.4f}')
col = [sum(r[1][i] for r in band_tot) for i in range(len(DEPTHS))]
print('  ' + '-' * (len(hdr) - 2))
print(f'  {"ALL":24s}' + ''.join(f'{v:+12.4f}' for v in col) + f'{sum(col):+10.4f}')

print('\n  ENTRY: mean net surface heat flux by band, W/m2 of GLOBAL area')
print('  (sign flipped from `fh` as stored, so positive = INTO the ocean;')
print('   this is entry, not storage, and circulation moves it afterwards)\n')
sa = area[0, :]
tot_e = 0.0
for name, la0, la1 in BANDS:
    sel = (lat >= la0) & (lat < la1)
    acc, n = 0.0, 0
    for y in LATE:
        f = glob.glob(f'{root}/outdata/fesom/fh.fesom.{y}.nc')
        if not f:
            continue
        with xr.open_dataset(f[0], decode_times=False) as d:
            nm = 'fh' if 'fh' in d.data_vars else list(d.data_vars)[-1]
            v = d[nm].mean(dim='time').values
        acc += float(np.nansum(-v[sel] * sa[sel])); n += 1   # negate: see SIGN note
    e = acc / n / A_EARTH if n else np.nan
    tot_e += e if np.isfinite(e) else 0.0
    print(f'  {name:24s} {e:+10.4f}')
print(f'  {"ALL":24s} {tot_e:+10.4f}')
