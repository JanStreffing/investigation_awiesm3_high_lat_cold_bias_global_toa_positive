"""Does the IFS sea-ice column reach its melting cap in summer?

With ECE_CPL_FESOM_ICE_QOUT the heat FESOM receives is base conduction plus the heat the
IFS column loses at its melting cap (RTMELTSICE = 273.15 K).  If the column's top layer
(istl1) and skin (skt) sit below melting through summer, almost no melt energy reaches
FESOM, whatever the coupling.  Cells: model ci >= 0.90, ocean, 60-90N, cos-lat weighted.
OIFS monthly remapped output.

Usage:  RUNS=11Y:1350-1350,15A:1350-1350 python3 scripts/analysis/ifs_ice_column_summer.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
LSMF=('/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc')
RUNS=os.environ.get('RUNS','11Y:1350-1350,15A:1350-1350,13A:1357-1359')
MN=['Apr','May','Jun','Jul','Aug','Sep']; MI=[3,4,5,6,7,8]
with xr.open_dataset(LSMF,decode_times=False) as d:
    lsm=np.squeeze(d['lsm'].values); lsm=lsm[0] if lsm.ndim==3 else lsm
    lat=np.squeeze(d['lat'].values)
def mload(arm,v,y):
    for p in (f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc',
              f'{R}/{arm}/run_13500101-13591231/work/atm_remapped_1m_{v}_{y}-{y}.nc'):
        if os.path.exists(p):
            with xr.open_dataset(p,decode_times=False) as d:
                k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
                return np.squeeze(d[k].values)
L2=np.broadcast_to(lat[:,None],lsm.shape); W=np.cos(np.deg2rad(L2)); base=(L2>=60)&(lsm<=0.5)
for item in RUNS.split(','):
    arm,yr=item.split(':'); y0,y1=map(int,yr.split('-'))
    acc={k:[[] for _ in MI] for k in ('istl1','skt','fmelt')}
    for y in range(y0,y1+1):
        ci,ist,skt=mload(arm,'ci',y),mload(arm,'istl1',y),mload(arm,'skt',y)
        if ci is None or ist is None or skt is None: continue
        for j,m in enumerate(MI):
            k=base&(ci[m]>=0.9)&np.isfinite(ist[m])
            if not k.any(): continue
            acc['istl1'][j].append(np.average(ist[m][k],weights=W[k]))
            acc['skt'][j].append(np.average(skt[m][k],weights=W[k]))
            acc['fmelt'][j].append(np.average((ist[m][k]>=273.0).astype(float),weights=W[k]))
    f=lambda L:[np.mean(x) if x else np.nan for x in L]
    print(f'\n{arm} {yr}   consolidated NH ice (ci>=0.9)   '+' '.join(f'{m:>7}' for m in MN))
    print('  istl1 K  (IFS top ice layer, cap 273.15) '+' '.join(f'{x:7.2f}' for x in f(acc['istl1'])))
    print('  skt K    (skin)                          '+' '.join(f'{x:7.2f}' for x in f(acc['skt'])))
    print('  cells with istl1 >= 273.0 (monthly mean) '+' '.join(f'{100*x:6.1f}%' for x in f(acc['fmelt'])))
