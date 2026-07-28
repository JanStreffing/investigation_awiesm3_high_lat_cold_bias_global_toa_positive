# Recovered from the Claude Code session transcript on 2026-07-28.
# Original location was the ephemeral session scratchpad, which was wiped.
# Session c56eada6-56b5-4b15-83cd-c6ed69cf48d9, written 07-27T14:24:58.

import numpy as np, xarray as xr, glob
files = sorted(glob.glob('/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_forcing/outdata/oifs/atm_1d_18*.nc'))
lsm = None
try:
    import cfgrib
    g = xr.open_dataset('/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_forcing/input/oifs/95_4/lsm', engine='cfgrib')
    lsm = g['lsm'].values
except Exception as e:
    print('no lsm:', e)

boxes = {'Boreal 55-70N':(55,70,-180,180), 'Siberia':(55,75,60,180), 'E.Siberia':(55,75,90,160)}
res={}
for f in files:
    ds = xr.open_dataset(f)
    lat = ds['lat'].values; lon = ((ds['lon'].values+180)%360)-180
    tas = ds['tas']
    mon = tas['time_counter.month'].values
    for name,(la0,la1,lo0,lo1) in boxes.items():
        m = (lat>=la0)&(lat<=la1)&(lon>=lo0)&(lon<=lo1)
        if lsm is not None: m &= (lsm>0.5)
        w = np.cos(np.deg2rad(lat[m]))
        jja = tas.values[np.isin(mon,[6,7,8])][:,m]
        ann = tas.values[:,m]
        res.setdefault(name,[]).append((np.average(jja.mean(0),weights=w)-273.15,
                                        np.average(ann.mean(0),weights=w)-273.15))
    ds.close()
print(f"{'box':16s} {'JJA mean':>9s} {'JJA sd':>7s} {'SE(8yr)':>8s} {'ANN mean':>9s} {'ANN sd':>7s}")
for name,v in res.items():
    a=np.array(v); jja=a[:,0]; ann=a[:,1]
    print(f"{name:16s} {jja.mean():9.2f} {jja.std(ddof=1):7.2f} {jja.std(ddof=1)/np.sqrt(8):8.2f} {ann.mean():9.2f} {ann.std(ddof=1):7.2f}")
    print(f"{'':16s} per-year JJA: {np.round(jja,2)}")
