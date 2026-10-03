"""Global-mean sea surface height over time, single vs double precision FESOM.

Under linfs the tracer volume is fixed, and the vertical velocity at the free surface
carries tracer through the surface (oce_adv_tra_ver.F90: tvert(nzmin) = -W*T*area).  If
single precision biases that velocity relative to the actual free-surface tendency, the
result is a spurious surface heat flux -- a candidate for the ~0.5 W/m2 SP ocean heat loss.
This measures whether the global volume (area-weighted mean ssh) drifts in SP and not DP.

Usage:  python3 scripts/analysis/ssh_volume_drift_sp_dp.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'
with xr.open_dataset('/work/ab0246/a270092/input/fesom2/core3/mesh.nc',decode_times=False) as m:
    A=m['cell_area'].values.astype('f8')
for arm,yrs in (('11X',range(1350,1390)),('11Y',range(1350,1390)),('13A',range(1350,1360))):
    v=[]
    for y in yrs:
        p=f'{R}/{arm}/outdata/fesom/ssh.fesom.{y}.nc'
        if not os.path.exists(p): break
        with xr.open_dataset(p,decode_times=False) as d:
            s=np.asarray(d['ssh'].values,dtype='f8').mean(0)
        k=np.isfinite(s); v.append(float((s[k]*A[k]).sum()/A[k].sum()))
    v=np.array(v)
    if len(v)<2: print(arm,'missing'); continue
    tr=np.polyfit(np.arange(len(v)),v,1)[0]
    print(f'{arm}: {len(v)} yr, global-mean ssh first {v[0]*100:+.2f} cm, last {v[-1]*100:+.2f} cm, '
          f'trend {tr*100:+.3f} cm/yr   (0.5 W/m2 via w_top*T would need ~30 cm/yr)')
