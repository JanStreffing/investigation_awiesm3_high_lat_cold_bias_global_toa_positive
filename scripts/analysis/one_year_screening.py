"""Can a Southern Ocean cloud lever be screened on ONE year instead of 44?

WHY THIS MATTERS FOR COST.  Project issue #170 found that the cy48 high-cloud offset is
already present in year 1 while the T2m bias and the SO low-cloud dip are not, and
concluded "one year long simulations will likely suffice to test the impact of the tuning
parameters".  A 1-year AMIP leg is ~11 minutes; the campaign's 44-year standard is ~9.5
hours.  If cloud-radiative response screens at 1 year, candidate levers can be triaged
50x cheaper and only survivors run to full length.

BUT THE CLAIM IS ABOUT hcc, AND MY METRIC IS SO SW CRE.  Those are different quantities
with different noise, so the transfer has to be measured, not assumed.  This script asks
directly: for levers whose 44-year answer is known, does year 1 alone give the same sign
and a comparable magnitude, and is the response bigger than the 1-year noise?

METHOD.  For each lever, SO 65-45S SW CRE delta against the shared control, computed (a)
from year 1 only and (b) from all 44 years.  The 1-year detection threshold comes from the
interannual scatter of the control -- the same discipline the campaign applies to
temperature, applied here before trusting any short run.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[_v]='1'
import numpy as np, xarray as xr, warnings; warnings.filterwarnings('ignore')
RT='/work/bb1469/a270092/runtime/oifsamip-cy48'; ACC=3600.0; SO=(-65.,-45.)
Y0,Y1=1872,1915
LEV=[('control','amip_pi_base'),('A1a ovlliqice0.10','amip_A1_overlap01'),
     ('A1b ovlliqice0.35','amip_A1_overlap035'),('D2a inpsea0.2','amip_D2a_inpsea02'),
     ('D2b inp+p700','amip_D2b_inpsea02_p700'),('B3 clddiff1.5e-5','amip_B3_clddiff15e6'),
     ('B7 rvice0.22','amip_B7_rvice022')]

def yearly_socre(run):
    out={}
    for y in range(Y0,Y1+1):
        fs=f'{RT}/{run}/outdata/oifs/atm_remapped_1m_tsr_1m_{y}-{y}.nc'
        fc=f'{RT}/{run}/outdata/oifs/atm_remapped_1m_tsrc_1m_{y}-{y}.nc'
        if not (os.path.exists(fs) and os.path.exists(fc)): continue
        with xr.open_dataset(fs,decode_times=False) as d:
            a=d['tsr'].values.mean(axis=0)/ACC; lat=d['lat'].values
        with xr.open_dataset(fc,decode_times=False) as d:
            b=d['tsrc'].values.mean(axis=0)/ACC
        s=(lat>=SO[0])&(lat<SO[1]); w=np.cos(np.deg2rad(lat[s]))
        out[y]=float(np.average((a-b)[s,:].mean(axis=1),weights=w))
    return out

print(__doc__); print('='*94)
ctl=yearly_socre('amip_pi_base')
yrs=sorted(ctl)
sd=np.std([ctl[y] for y in yrs],ddof=1)
thr1=1.96*sd*np.sqrt(2.0)          # 95 % on a difference of two SINGLE years
thr44=1.96*sd*np.sqrt(2.0/len(yrs))
print(f'control SO SW CRE: {len(yrs)} yr, mean {np.mean([ctl[y] for y in yrs]):.2f}, '
      f'interannual sd {sd:.3f} W/m2')
print(f'  95 % detection threshold:  1-year pair +-{thr1:.2f}   44-year pair +-{thr44:.2f}\n')
print(f'{"lever":20s} {"yr1 delta":>10s} {"44yr delta":>11s} {"ratio":>7s}  '
      f'{"yr1 signif":>10s}  verdict')
for nm,run in LEV[1:]:
    d=yearly_socre(run)
    common=sorted(set(d)&set(ctl))
    if not common: continue
    y1=d[yrs[0]]-ctl[yrs[0]] if yrs[0] in d else np.nan
    full=np.mean([d[y]-ctl[y] for y in common])
    sig='YES' if abs(y1)>thr1 else 'no'
    agree = (np.sign(y1)==np.sign(full)) and abs(y1)>thr1
    print(f'{nm:20s} {y1:+10.2f} {full:+11.2f} {y1/full if full else np.nan:7.2f}  '
          f'{sig:>10s}  {"screens at 1 yr" if agree else "NOT reliable at 1 yr"}')
print(f'\n  A lever screens at 1 year only if its year-1 delta clears +-{thr1:.2f} AND has')
print('  the same sign as the 44-year answer.  Anything smaller needs the full length.')
