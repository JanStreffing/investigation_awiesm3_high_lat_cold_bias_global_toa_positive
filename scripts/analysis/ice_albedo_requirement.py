"""How much brighter must the sea ice be, and on WHICH namelist knob?

Two questions, kept separate.

1. HOW MUCH ENERGY.  To carry dV more ice through to the summer minimum, the melt season
   must absorb dV * rho_ice * L_f fewer joules.  Spread over the melt-season pack area and
   length that is a W/m2, and divided by the pack's downwelling shortwave it is an albedo
   change.  This is a LOWER BOUND on the albedo needed: it credits no feedback, and the
   ice-albedo feedback amplifies whatever is applied.

2. WHICH KNOB.  FESOM's ice_albedo (ice_thermo_cpl.F90:701-738) is, per ice node,

       alb_snow = albsn + fmelt*(albsnm - albsn)      fmelt = clamp((T-273.15+tramp)/tramp)
       alb_bare = albi  + fmelt*(albim  - albi)
       alb_noponds = alb_snow if hsn > 1 mm else alb_bare        (h_snowscale = 0 -> STEP)
       alb = (1-apnd)*alb_noponds + apnd*albpnd                  (meltpond_albedo, aice=1)

   so the sensitivity of the pack albedo to each of the four knobs is an area fraction,
   measurable from the run.  d(alb)/d(albsnm) = mean[(1-apnd) * 1{hsn>1mm} * fmelt], and so
   on.  h_snowscale = 0 means the snow/bare split is a HARD STEP at 1 mm, which is the
   nonlinearity behind the 2026-09 runaway: brightening melting snow stops the pack
   clearing its snow, which keeps it on the snow branch, which brightens it further.
   The NH snow-free fraction is printed because it is the margin against that.

Usage:  ARM=PICAL_ccnice Y0=1960 Y1=1969 python3 scripts/analysis/ice_albedo_requirement.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R   = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PICAL_ccnice')
Y0, Y1 = int(os.environ.get('Y0', 1960)), int(os.environ.get('Y1', 1969))
RHO, LF = 917.0, 3.34e5          # kg/m3, J/kg
TRAMP = 1.0                      # alb_tramp, staged namelist.ice
FG = f'{R}/{ARM}/outdata/fesom'
OG = f'{R}/{ARM}/outdata/oifs'

def gr(var, y):
    p = f'{FG}/{var}.fesom.gr.{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p, decode_times=False) as d:
        return np.asarray(d[var].values, np.float32), np.asarray(d['lat'].values, float)

def months_of(n):
    """day index -> month index; handles 12, 365 and leap 366."""
    if n == 12: return np.arange(12)
    ml = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    if n == 366: ml = ml.copy(); ml[1] = 29
    return np.repeat(np.arange(12), ml)[:n]

# melt seasons: NH Jun-Aug, SH Dec-Feb
SEAS = {'NH': ([5,6,7],  lambda la: la > 45), 'SH': ([11,0,1], lambda la: la < -45)}

print(__doc__.split('\n')[0]); print(f'{ARM}, melt seasons of {Y0}-{Y1}\n')

acc = {h: {k: [] for k in ('d_albsnm','d_albim','d_albsn','d_albi',
                           'snowfree','pond','area','vol')} for h in SEAS}
for y in range(Y0, Y1+1):
    A = gr('a_ice', y); S = gr('m_snow', y); P = gr('apnd', y); T = gr('ist', y)
    if A is None or S is None or P is None or T is None: continue
    a, lat = A; sn, _ = S; pnd, _ = P; ist, _ = T
    # apnd/ist are NaN off-ice; the pack mask handles that, but the weighted means
    # need finite values inside it too.
    pnd = np.nan_to_num(pnd); a = np.nan_to_num(a); sn = np.nan_to_num(sn)
    ist = ist + 273.15 if np.nanmean(ist) < 100 else ist
    dlat = abs(lat[1]-lat[0])
    cell = (6.371e6**2 * np.cos(np.deg2rad(lat)) * np.deg2rad(dlat)
            * np.deg2rad(360.0/a.shape[2]))[:, None]
    ma, mp = months_of(a.shape[0]), months_of(pnd.shape[0])
    mt = months_of(ist.shape[0])
    for h, (mons, latsel) in SEAS.items():
        ka = np.isin(ma, mons); kp = np.isin(mp, mons); kt = np.isin(mt, mons)
        aa = np.nanmean(a[ka], 0); ss = np.nanmean(sn[ka], 0)
        pp = np.nanmean(pnd[kp], 0); tt = np.nanmean(ist[kt], 0)
        band = np.broadcast_to(latsel(lat)[:, None], aa.shape)
        pack = band & (aa > 0.15) & np.isfinite(aa)
        if not pack.any(): continue
        w = np.broadcast_to(cell, aa.shape)[pack] * aa[pack]     # ice-area weighted
        hsn = np.where(aa > 0.01, ss/np.maximum(aa, 0.01), 0.0)  # snow depth ON the ice
        snow = (hsn[pack] > 0.001)
        fm   = np.clip((np.nan_to_num(tt[pack], nan=273.15) - 273.15 + TRAMP)/TRAMP, 0, 1)
        op   = 1.0 - np.clip(pp[pack], 0, 1)
        acc[h]['d_albsnm'].append(np.average(op*snow*fm, weights=w))
        acc[h]['d_albim' ].append(np.average(op*(~snow)*fm, weights=w))
        acc[h]['d_albsn' ].append(np.average(op*snow*(1-fm), weights=w))
        acc[h]['d_albi'  ].append(np.average(op*(~snow)*(1-fm), weights=w))
        acc[h]['snowfree'].append(1.0 - np.average(snow, weights=w))
        acc[h]['pond'    ].append(np.average(np.clip(pp[pack],0,1), weights=w))
        acc[h]['area'    ].append(float((np.broadcast_to(cell, aa.shape)[pack]*aa[pack]).sum()))

# ---- melt-season downwelling SW over the pack, from OIFS (accumulated J/m2 -> W/m2)
def oifs(var, y):
    p = f'{OG}/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(p): return None
    with xr.open_dataset(p, decode_times=False) as d:
        k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        return np.squeeze(d[k].values), np.squeeze(d['lat'].values)
SW = {}
for h, (mons, latsel) in SEAS.items():
    v = []
    for y in range(Y0, Y1+1):
        D = oifs('ssrd', y); C = oifs('ci', y)
        if D is None or C is None: continue
        d_, la = D; c_, _ = C
        band = np.broadcast_to(latsel(la)[:, None], d_.shape[1:])
        k = band & (c_[mons].mean(0) > 0.15)
        if k.any(): v.append(float(np.nanmean(d_[mons].mean(0)[k]))/3600.0)
    SW[h] = float(np.mean(v)) if v else np.nan

print(f'{"":<6}{"SWdn":>8}{"snow-free":>11}{"pond":>8}'
      f'{"d/dalbsnm":>11}{"d/dalbim":>10}{"d/dalbsn":>10}{"d/dalbi":>9}{"pack Mkm2":>11}')
SENS = {}
for h in SEAS:
    m = {k: float(np.mean(v)) if v else np.nan for k, v in acc[h].items()}
    SENS[h] = m; SENS[h]['sw'] = SW[h]
    print(f'{h:<6}{SW[h]:>8.1f}{m["snowfree"]*100:>10.1f}%{m["pond"]*100:>7.1f}%'
          f'{m["d_albsnm"]:>11.3f}{m["d_albim"]:>10.3f}{m["d_albsn"]:>10.3f}'
          f'{m["d_albi"]:>9.3f}{m["area"]/1e12:>11.2f}')

print('\nenergy requirement (lower bound, no feedback credit):')
print(f'{"":<6}{"target dV":>12}{"dE":>12}{"dF":>10}{"d(alb) needed":>15}')
TARG = {'NH': 4.0e12, 'SH': 2.0e12}     # m3 to recover at the summer minimum
TMELT = 92*86400.0
for h in SEAS:
    dE = TARG[h]*RHO*LF
    dF = dE/(SENS[h]['area']*TMELT)
    da = dF/SENS[h]['sw']
    SENS[h]['da'] = da
    print(f'{h:<6}{TARG[h]/1e12:>10.1f}e3km3{dE:>12.2e}{dF:>10.2f}{da:>15.3f}')
print('\n  dF is W/m2 of melt-season shortwave to be reflected instead of absorbed.')
print('  Feedback amplifies, so the applied change should be BELOW these numbers.')
