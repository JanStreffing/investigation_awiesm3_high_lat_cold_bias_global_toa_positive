"""Net TOA next to the snow-enthalpy term, per year, for one coupled arm.

The snow enthalpy (latent heat of fusion carried by snowfall, sf * rho_w * Lf =
sf * 3.3355e8 J/m of water) is not part of TOA.  It belongs in the SURFACE budget: snow
that falls into the ocean is melted by the ocean, so an ocean heat uptake measured against
TOA or against the IFS surface flux differs from it by this term.  Printed as a global mean
and as the part that falls on ocean, both in W/m2 of Earth area.

Usage:  ARM=15F Y0=1360 Y1=1369 python3 scripts/analysis/toa_vs_snow_enthalpy.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'; ACC=3600.0; LF_RHO=3.3355e8
ARM=os.environ.get('ARM','15F'); Y0=int(os.environ.get('Y0',1360)); Y1=int(os.environ.get('Y1',1369))
LSMF=('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
with xr.open_dataset(LSMF,decode_times=False) as d:
    m=np.squeeze(d['lsm'].values); m=m[0] if m.ndim==3 else m
    lat=np.squeeze(d['lat'].values)
W=np.cos(np.deg2rad(np.broadcast_to(lat[:,None],m.shape))); ocean=(m<=0.5)
def ld(v,y):
    with xr.open_dataset(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc',decode_times=False) as d:
        k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]; return np.squeeze(d[k].values)
gm=lambda f: float((f*W).sum()/W.sum())          # global mean, W/m2 of Earth
print(f'{ARM} {Y0}-{Y1}   [W/m2 of Earth area]')
print(f'  {"year":>5}{"net TOA":>10}{"snow enth":>11}{"  on ocean":>11}{"TOA - ocean snow":>18}')
rows=[]
for y in range(Y0,Y1+1):
    toa=gm((ld('tsr',y)+ld('ttr',y)).mean(0)/ACC)
    sf=ld('sf',y).mean(0)/ACC*LF_RHO                 # accumulated m w.e. per hour -> W/m2
    se=gm(sf); so=gm(np.where(ocean,sf,0.0))
    rows.append((toa,se,so)); print(f'  {y:>5}{toa:10.3f}{se:11.3f}{so:11.3f}{toa-so:18.3f}')
a=np.array(rows).mean(0)
print(f'  {"mean":>5}{a[0]:10.3f}{a[1]:11.3f}{a[2]:11.3f}{a[0]-a[2]:18.3f}')
