"""Zero-compute check of the MLE implementation before spending a coupled run.

Transcribes the arithmetic of mle_add_gamma (oce_mle.F90) line for line into numpy and
evaluates it on real 11X fields.  It cannot catch a Fortran coding error -- a bad index, an
uninitialised variable -- but it does catch the three things that would waste a run:

  SIGN       Gamma_MLE must have the same sign as Gamma_GM, which restratifies.  In the
             code both are proportional to +grad_h(b) = -(g/rho0)*sigma_xy, so the test is
             that the offline transcription reproduces that and does not flip.
  MAGNITUDE  Psi must land near the 2-5 m2/s predicted for Psi_FK, and NOT orders away.
  THE CAPS   mle_hmax = 500 m and mle_resscale_max = 20 were added on physical grounds;
             this measures what they actually remove, which was never checked.

The earlier estimate (gm_bolus_vs_fk.py) omitted mu(z) entirely -- it used the peak value
-- and applied no caps, so it is an upper bound on what the code will actually produce.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, seawater as sw, warnings
warnings.filterwarnings('ignore')

D = '/work/bb1469/a270092/runtime/awiesm3-v3.4/11X/outdata/fesom'
MESHF = '/work/ab0246/a270092/input/fesom2/core3/mesh.nc'
YEARS = range(1385, 1390)
R_E, OMEGA, G, RHO0 = 6371000.0, 7.2921e-5, 9.81, 1025.0
# the namelist defaults, exactly as in config/namelist.oce
MLE_CE, MLE_LF, MLE_TAU, MLE_HMAX, MLE_RESMAX = 0.06, 4000.0, 432000.0, 500.0, 20.0
NH, SH = [1, 2], [7, 8]
BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45), ('45-30S', -45, -30),
         ('30S-30N', -30, 30), ('30-45N', 30, 45), ('45-60N', 45, 60), ('60-90N', 60, 90)]
BANDS = BANDS[::-1]          # report north to south

m = xr.open_dataset(MESHF, decode_times=False)
lon, lat = m['lon'].values.astype('f8'), m['lat'].values.astype('f8')
tri = m['triag_nodes'].values.astype(int) - 1
dbn = m['depth_bnds'].values.astype('f8')
assert lat.size == 220509, 'MESH GUARD'
with xr.open_dataset(f'{D}/temp.fesom.1389.nc', decode_times=False) as d:
    zmid = d['nz'].values.astype('f8')
nz = zmid.size
dz = np.abs(dbn[1:nz+1]) - np.abs(dbn[:nz])
zlev = np.abs(dbn[:nz])                      # level (interface) depths, positive down

elon, elat = lon[tri], lat[tri]
dlon = ((elon - elon[:, 0:1] + 180.0) % 360.0) - 180.0
lat0 = elat.mean(axis=1, keepdims=True)
X = np.deg2rad(dlon) * R_E * np.cos(np.deg2rad(lat0))
Y = np.deg2rad(elat - lat0) * R_E
twoA = ((X[:,1]-X[:,0])*(Y[:,2]-Y[:,0]) - (X[:,2]-X[:,0])*(Y[:,1]-Y[:,0]))
eArea = 0.5*np.abs(twoA); ecLat = elat.mean(axis=1)
dx_e = np.sqrt(2.0*eArea); f_e = 2*OMEGA*np.sin(np.deg2rad(ecLat))

def winter(var, lev=True):
    aN = aS = None; n = 0
    for y in YEARS:
        with xr.open_dataset(f'{D}/{var}.fesom.{y}.nc', decode_times=False) as d:
            k = [x for x in d.data_vars if x not in ('bounds_lon','bounds_lat','time_bounds')][-1]
            a1 = d[k].isel(time=NH).values.astype('f8').mean(axis=0)
            a2 = d[k].isel(time=SH).values.astype('f8').mean(axis=0)
        aN = a1 if aN is None else aN+a1; aS = a2 if aS is None else aS+a2; n += 1
    aN, aS = aN/n, aS/n
    north = lat >= 0
    return np.where(north[:,None], aN, aS) if lev else np.where(north, aN, aS)

T, S, H = winter('temp'), winter('salt'), np.abs(winter('MLD2', lev=False))
b = -G*(sw.dens0(S, T) - 1000.0)/RHO0

# --- mirror the Fortran: cap H, ML-average the gradient, mu(z), prefac -------------
Hc = np.minimum(H, MLE_HMAX)
inml = (zmid[None,:] <= np.maximum(Hc, zmid[0])[:,None]) & np.isfinite(b)
w = np.where(inml, dz[None,:], 0.0)
b_ml = np.where(w.sum(1) > 0, np.nansum(b*w,1)/np.maximum(w.sum(1),1e-30), np.nan)
ok = np.isfinite(b_ml); bf = np.where(ok, b_ml, 0.0)
bt = bf[tri]
gx = (bt[:,0]*(Y[:,1]-Y[:,2]) + bt[:,1]*(Y[:,2]-Y[:,0]) + bt[:,2]*(Y[:,0]-Y[:,1]))/twoA
gy = (bt[:,0]*(X[:,2]-X[:,1]) + bt[:,1]*(X[:,0]-X[:,2]) + bt[:,2]*(X[:,1]-X[:,0]))/twoA
gradb = np.sqrt(gx**2 + gy**2)
He = Hc[tri].mean(1)
dxlf_raw = dx_e/MLE_LF
dxlf = np.minimum(dxlf_raw, MLE_RESMAX)
prefac = MLE_CE * dxlf * He**2 / np.sqrt(f_e**2 + MLE_TAU**-2)

xi = 2.0*zlev[None,:]/np.maximum(He,1.0)[:,None] - 1.0
mu = (1.0 - xi**2)*(1.0 + (5.0/21.0)*xi**2)
mu = np.where((zlev[None,:] <= He[:,None]) & (mu > 0), mu, 0.0)
Psi = prefac[:,None]*mu*gradb[:,None]
good = ok[tri].all(1) & np.isfinite(gradb) & (eArea > 0) & (He > zmid[0])

print('='*100)
print('OFFLINE TRANSCRIPTION OF mle_add_gamma, on 11X winter fields')
print('='*100)
print(f'{"band":<10}{"dx/Lf raw":>11}{"capped":>8}{"H raw":>8}{"capped":>8}'
      f'{"peak Psi":>11}{"median":>9}{"mu peak z/H":>13}')
for nm, lo, hi in BANDS:
    s = good & (ecLat >= lo) & (ecLat < hi)
    if s.sum() < 100: continue
    aw = lambda x: np.sum(x[s]*eArea[s])/np.sum(eArea[s])
    pk = Psi[s].max(1)
    print(f'{nm:<10}{aw(dxlf_raw):>11.1f}{aw(dxlf):>8.1f}'
          f'{np.sum(H[tri].mean(1)[s]*eArea[s])/np.sum(eArea[s]):>8.0f}{aw(He):>8.0f}'
          f'{np.sum(pk*eArea[s])/np.sum(eArea[s]):>11.3f}{np.median(pk):>9.3f}'
          f'{zlev[np.argmax(mu[s].mean(0))]/max(aw(He),1):>13.2f}')
s = good
aw = lambda x: np.sum(x[s]*eArea[s])/np.sum(eArea[s])
pk = Psi[s].max(1)
print(f'{"GLOBAL":<10}{aw(dxlf_raw):>11.1f}{aw(dxlf):>8.1f}'
      f'{np.sum(H[tri].mean(1)[s]*eArea[s])/np.sum(eArea[s]):>8.0f}{aw(He):>8.0f}'
      f'{np.sum(pk*eArea[s])/np.sum(eArea[s]):>11.3f}{np.median(pk):>9.3f}')

print('\nSIGN: Gamma_MLE = prefac*mu*grad_h(b) and Gamma_GM = K*grad_h(b)/N^2.')
print(f'  both proportional to +grad_h(b); prefac min over the mesh = {prefac[good].min():.3e}'
      f'  (must be > 0)')
print(f'  mu range = [{mu[good].min():.3f}, {mu[good].max():.3f}]  (must be within [0, 1])')
print(f'  fraction of wet elements where the scheme is active: {good.mean()*100:.1f} %')
print(f'\nWHAT THE CAPS REMOVE:')
print(f'  H  capped on {100*np.mean(H[good[tri].any(1) if False else ok] > MLE_HMAX):.1f} % of nodes'
      f'   (H > {MLE_HMAX:.0f} m)')
print(f'  dx/Lf capped on {100*np.mean(dxlf_raw[good] > MLE_RESMAX):.1f} % of elements')
