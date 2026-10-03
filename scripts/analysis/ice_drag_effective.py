"""Effective atmosphere-ice drag in the coupled model, SH winter (JJA), per sector.
Cd_ai = |<tau_ai>| / (rho_a |<|U10| U10>|): monthly-mean FESOM atmice stress (what OpenIFS hands over)
against the monthly mean of the daily vector |U|U from OpenIFS 10 m wind, so a quadratic drag law with a
constant Cd is recovered exactly. Pack only (monthly a_ice >= 0.8), 55-80S. Also the free-drift wind factor
implied by Cd_ai and FESOM's Cd_oce_ice: sqrt(rho_a Cd_ai / (rho_w Cd_io)), and mean ice strength and mass.
Usage: python ice_drag_effective.py <run outdata> <y0> <y1> <Cd_oce_ice>
"""
import sys, numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R, y0, y1, cdio = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
RA, RW = 1.3, 1027.0
acc = {}
for y in range(y0, y1 + 1):
    f = lambda v: xr.open_dataset(f'{R}/fesom/{v}.fesom.gr.{y}.nc', decode_times=False)[v]
    tx, ty, a, st, mi = f('atmice_x'), f('atmice_y'), f('a_ice'), f('strength_ice'), f('m_ice')
    a = a if a.sizes['time'] == 12 else None
    U = xr.open_dataset(f'{R}/oifs/atm_remapped_1d_10u_{y}-{y}.nc', use_cftime=True)['10u']
    V = xr.open_dataset(f'{R}/oifs/atm_remapped_1d_10v_{y}-{y}.nc', use_cftime=True)['10v']
    td = 'time_counter' if 'time_counter' in U.dims else 'time'
    sp = np.hypot(U, V); su, sv = (sp * U).groupby(f'{td}.month').mean(), (sp * V).groupby(f'{td}.month').mean()
    ad = xr.open_dataset(f'{R}/fesom/a_ice.fesom.gr.{y}.nc', decode_times=False)['a_ice']
    for m in (6, 7, 8):
        ai = ad.isel(time=slice(sum([31,28,31,30,31,30,31,31][:m-1]), sum([31,28,31,30,31,30,31,31][:m]))).mean('time') if ad.sizes['time'] > 12 else ad.isel(time=m-1)
        suu = su.sel(month=m).interp(lat=tx.lat, lon=tx.lon % 360 if float(su.lon.min()) >= 0 else tx.lon)
        svv = sv.sel(month=m).interp(lat=tx.lat, lon=tx.lon % 360 if float(su.lon.min()) >= 0 else tx.lon)
        rec = dict(tau=np.hypot(tx.isel(time=m-1), ty.isel(time=m-1)).values, uu=np.hypot(suu, svv).values,
                   ai=ai.values, st=st.isel(time=m-1).values, mi=mi.isel(time=m-1).values)
        for k, v in rec.items(): acc.setdefault(k, []).append(v)
A = {k: np.stack(v) for k, v in acc.items()}
lat, lon = tx.lat.values, ((tx.lon.values + 180) % 360) - 180
LON, LAT = np.meshgrid(lon, lat); w = np.cos(np.deg2rad(LAT))
SECT = {'all': np.ones_like(LON, bool), 'Weddell': (LON >= -60) & (LON < 20), 'Indian': (LON >= 20) & (LON < 90),
        'Pacific': (LON >= 90) & (LON < 160), 'Ross': (LON >= 160) | (LON < -130), 'Amund-Bell': (LON >= -130) & (LON < -60)}
print(f'SH JJA {y0}-{y1}, pack a_ice >= 0.8, 55-80S.  Cd_oce_ice = {cdio}')
print(f'{"sector":11s} {"|tau| N/m2":>10s} {"|<|U|U>| m2/s2":>14s} {"Cd_ai x1e3":>10s} {"free-drift factor %":>19s} {"strength N/m":>12s} {"m_ice m":>8s}')
for s, sel in SECT.items():
    m = sel[None] & (A['ai'] >= 0.8) & (LAT[None] < -55) & (LAT[None] > -80) & np.isfinite(A['tau']) & np.isfinite(A['uu']) & (A['uu'] > 1)
    W = np.broadcast_to(w, m.shape)[m]
    tau = np.average(A['tau'][m], weights=W); uu = np.average(A['uu'][m], weights=W)
    cd = tau / (RA * uu)
    print(f'{s:11s} {tau:10.3f} {uu:14.1f} {cd*1e3:10.2f} {100*np.sqrt(RA*cd/(RW*cdio)):19.2f} {np.average(A["st"][m], weights=W):12.0f} {np.average(A["mi"][m], weights=W):8.2f}')
