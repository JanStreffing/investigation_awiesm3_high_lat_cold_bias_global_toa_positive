# Recovered from the Claude Code session transcript on 2026-07-28.
# Original location was the ephemeral session scratchpad, which was wiped.
# Session c56eada6-56b5-4b15-83cd-c6ed69cf48d9, written 07-28T08:24:30.

import numpy as np, xarray as xr, os, warnings; warnings.filterwarnings('ignore')
ACC=3600.0; SFCONV=333_550_000.0
def gm(ds,v):
    da=ds[v].mean('time_counter')/ACC
    w=np.cos(np.deg2rad(da.lat.values))[:,None]*np.ones(da.shape)
    return float((da.values*w).sum()/w.sum())
def run(tag,d,inf,y0,y1):
    acc={k:[] for k in ['tsr','ttr','ssr','str','sshf','slhf','sf']}
    for y in range(y0,y1+1):
        ok=True; vals={}
        for v in acc:
            f=f'{d}/atm_remapped_1m_{v}{inf}_{y}-{y}.nc'
            if not os.path.exists(f): ok=False; break
            ds=xr.open_dataset(f); vals[v]=gm(ds,v); ds.close()
        if ok:
            for v in acc: acc[v].append(vals[v])
    m={k:np.mean(v) for k,v in acc.items()}
    toa=m['tsr']+m['ttr']
    sfc_nosnow=m['ssr']+m['str']+m['sshf']+m['slhf']
    snow=m['sf']*SFCONV
    sfc=sfc_nosnow-snow
    print(f"{tag:16s} n={len(acc['tsr']):2d}  TOA={toa:7.3f}  SFC(no snow)={sfc_nosnow:7.3f}  "
          f"snow-enth={snow:6.3f}  SFC(part2)={sfc:7.3f}  TOA-SFC={toa-sfc:7.3f}")
A='/work/bb1469/a270092/runtime/oifsamip-cy48/amip_pi_base/outdata/oifs'
C='/work/bb1469/a270092/runtime/awiesm3-v3.4/Tuning_test_09C_06V_CRUNCEPinit_newSeaIce/outdata/oifs'
B='/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_080a_lpjguess_Baseline_coupled_fromCRUNCEP/outdata/oifs'
run('AMIP 1872-79',A,'_1m',1872,1879)
run('09C 1370-79',C,'',1370,1379)
run('080a 1370-79',B,'',1370,1379)
