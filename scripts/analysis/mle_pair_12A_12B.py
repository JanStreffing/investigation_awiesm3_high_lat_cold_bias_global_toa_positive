"""12A (Fox-Kemper MLE on) against 12B (off): does it restratify, and by how much?

THE PAIR.  One year each from the same core3_linfs restart, differing in one namelist
entry.  12B ran library a4ac4d38 and 12A the limiter/level-range build; the MLE code is
inert when use_mle is false, so the arms are comparable, but that is an assumption worth
stating rather than hiding.

THE ONE QUESTION THIS RUN ANSWERS.  Is the SIGN right?  The mixed layer must SHOAL.  If it
deepens, Gamma_MLE is inverted and no amount of tuning helps -- stop and fix.

Magnitude is secondary here and one year cannot score it, but the added streamfunction is
checked against what was predicted offline:

    d(Psi) = |int bolus dz to the ML base|_12A - |...|_12B

isolates the MLE part because GM is common to both arms.  Predicted after the u* limiter:
0.9-2.0 m2/s at the outcrops, against a GM mixed-layer streamfunction of 0.10-0.23.

MLD is recomputed from T and S under WOA18's criterion (0.125 kg/m3 referenced to 10 m)
rather than read from MLD2, so the numbers are comparable with the baseline in
mld_baseline_vs_woa18.py.  MESH GUARD: core3, 220509 nodes.
"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, seawater as sw, warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
MESHF = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
YEAR = 1350
SIGMA_CRIT, REF_DEPTH = 0.125, 10.0
BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45), ('45-30S', -45, -30),
         ('30S-30N', -30, 30), ('30-45N', 30, 45), ('45-60N', 45, 60), ('60-90N', 60, 90)]
BANDS = BANDS[::-1]          # report north to south
MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

m = xr.open_dataset(MESHF, decode_times=False)
lat = m['lat'].values.astype('f8'); area = m['cell_area'].values.astype('f8')
tri = m['triag_nodes'].values.astype(int) - 1
dbn = m['depth_bnds'].values.astype('f8')
if lat.size != 220509: sys.exit('MESH GUARD: not core3')
with xr.open_dataset(f'{R}/12B/outdata/fesom/temp.fesom.{YEAR}.nc', decode_times=False) as d:
    zmid = d['nz'].values.astype('f8')
nz = zmid.size
dz = np.abs(dbn[1:nz+1]) - np.abs(dbn[:nz])
eArea = 0.5*np.abs((lat[tri][:,1]-lat[tri][:,0]))  # placeholder, replaced below
ecLat = lat[tri].mean(axis=1)

def get(arm, var):
    f = f'{R}/{arm}/outdata/fesom/{var}.fesom.{YEAR}.nc'
    if not os.path.exists(f): return None
    with xr.open_dataset(f, decode_times=False) as d:
        k = [x for x in d.data_vars if x not in ('bounds_lon','bounds_lat','time_bounds')][-1]
        return d[k].values.astype('f8')

def mld_woa(T, S):
    sig = sw.dens0(S, T) - 1000.0
    k = int(np.searchsorted(zmid, REF_DEPTH)); k = min(max(k,1), nz-1)
    w = (REF_DEPTH - zmid[k-1])/(zmid[k]-zmid[k-1])
    ref = sig[:,k-1] + w*(sig[:,k]-sig[:,k-1])
    ref = np.where(np.isfinite(ref), ref, sig[:,0])
    ex = sig - ref[:,None]; valid = ~np.isnan(sig)
    over = (ex >= SIGMA_CRIT) & valid
    out = np.full(sig.shape[0], np.nan); hit = over.any(1); first = np.argmax(over,1)
    idx = np.where(hit)[0]; j = first[idx]; jl = np.maximum(j-1,0)
    eh, el = ex[idx,j], ex[idx,jl]; den = eh-el
    fr = np.clip(np.where(np.abs(den)>1e-12,(SIGMA_CRIT-el)/den,0.0),0,1)
    out[idx] = zmid[jl] + fr*(zmid[j]-zmid[jl]); out[idx[j==0]] = zmid[0]
    nlev = valid.sum(1); wet = nlev>0
    bot = np.where(wet, zmid[np.clip(nlev-1,0,nz-1)], np.nan)
    out[wet & ~hit] = bot[wet & ~hit]; out[~wet] = np.nan
    return out

def band(fld, sel_lat, w):
    o=[]
    for nm,lo,hi in BANDS:
        s=(sel_lat>=lo)&(sel_lat<hi)&np.isfinite(fld)
        o.append(np.sum(fld[s]*w[s])/np.sum(w[s]) if s.any() else np.nan)
    return np.array(o)

TA,SA = get('12A','temp'), get('12A','salt')
TB,SB = get('12B','temp'), get('12B','salt')
if TA is None: sys.exit('12A has no output yet')

print('='*96)
print(f'MIXED LAYER DEPTH [m], WOA criterion, 12A (MLE on) minus 12B (off), year {YEAR}')
print('  NEGATIVE = the mixed layer SHOALED = the sign is right')
print('='*96)
print(f'{"month":<7}' + ''.join(f'{n:>11}' for n,_,_ in BANDS))
ann = np.zeros(len(BANDS)); nm_ = 0
for mth in range(TA.shape[0]):
    a = mld_woa(TA[mth], SA[mth]); b = mld_woa(TB[mth], SB[mth])
    d = band(a,lat,area) - band(b,lat,area)
    ann += d; nm_ += 1
    print(f'{MON[mth]:<7}' + ''.join(f'{v:>+11.2f}' for v in d))
print('-'*96)
print(f'{"MEAN":<7}' + ''.join(f'{v:>+11.2f}' for v in ann/max(nm_,1)))

# --- the added streamfunction -------------------------------------------------------
bua, bva = get('12A','bolus_u'), get('12A','bolus_v')
bub, bvb = get('12B','bolus_u'), get('12B','bolus_v')
if bua is not None and bub is not None:
    print('\n' + '='*96)
    print('ADDED BOLUS STREAMFUNCTION at the mixed layer base [m2/s]  (12A - 12B)')
    print('  predicted after the u* limiter: 0.9-2.0 at the outcrops; GM alone gives 0.10-0.23')
    print('='*96)
    eA = np.ones(tri.shape[0])
    print(f'{"month":<7}' + ''.join(f'{n:>11}' for n,_,_ in BANDS))
    for mth in [1,2,7,8]:
        Ha = mld_woa(TA[mth], SA[mth]); He = Ha[tri].mean(1)
        ca = np.nancumsum(np.nan_to_num(bua[mth])*dz[None,:],1)
        cb = np.nancumsum(np.nan_to_num(bub[mth])*dz[None,:],1)
        da = np.nancumsum(np.nan_to_num(bva[mth])*dz[None,:],1)
        db = np.nancumsum(np.nan_to_num(bvb[mth])*dz[None,:],1)
        k = np.clip(np.searchsorted(zmid, He), 0, nz-1); i = np.arange(len(k))
        pA = np.sqrt(ca[i,k]**2 + da[i,k]**2); pB = np.sqrt(cb[i,k]**2 + db[i,k]**2)
        print(f'{MON[mth]:<7}' + ''.join(f'{v:>+11.3f}' for v in band(pA-pB, ecLat, eA)))

mp = get('12A','mle_psi')
print(f"\nmle_psi in the output: {'yes' if mp is not None else 'no (def_stream under XIOS)'}")
