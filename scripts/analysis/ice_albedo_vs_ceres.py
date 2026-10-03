"""Surface albedo over the melt-season pack, model against CERES, and what a change would buy.

The question this answers: is the sea ice too dark, and by how much in W/m2 of absorbed
shortwave.  It is asked in the melt season of each hemisphere, because that is when the pack
sees any sun at all: NH June-August, SH December-February.

Cells are those where BOTH the model (ci) and HadISST (sic) carry >= 0.8 ice in that month, so
model and observation are compared over the same surface and an ice-extent error cannot
masquerade as an albedo error.  Albedo is the broadband surface value 1 - ssr/ssrd for the model
(IFS fluxes are accumulated J/m2 per output hour, but the ratio is insensitive to that) and
sfc_sw_up/sfc_sw_down for CERES.  `fal` is the albedo OpenIFS actually used, reported alongside
as a cross-check on the flux-derived number.

CAVEAT, and it is not a small one: CERES surface fluxes are a radiative-transfer product, not a
measurement, and they are at their weakest over ice - CERES and MODIS disagree by 7.7 pp on cloud
in 90-65S.  Read the model-CERES difference as a guide, and weigh the in-situ literature values
printed at the end at least as heavily.

Usage:  RUNS=PICAL_momixoff Y0=1920 Y1=1939 python3 scripts/analysis/ice_albedo_vs_ceres.py
"""
import os
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
RUNS = os.environ.get('RUNS', 'PICAL_momixoff').split(',')
Y0, Y1 = int(os.environ.get('Y0', 1920)), int(os.environ.get('Y1', 1939))
CER = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
SIC = '/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc'
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'
# CONC: the concentration both model and HadISST must exceed.  0.8 isolates solid pack, but in
# the SH melt season almost nothing reaches 0.8 in HadISST, so the southern answer needs a lower
# bar; run it at both and report both.
CONC = float(os.environ.get('CONC', 0.8))
SEASON = {'NH melt season (JJA)': ([5, 6, 7], 55, 90), 'SH melt season (DJF)': ([11, 0, 1], -90, -55),
          'SH shoulder (ON)': ([9, 10], -90, -55), 'NH shoulder (MJ)': ([4, 5], 55, 90)}

with xr.open_dataset(CER, decode_times=False) as d:
    clat = np.squeeze(d['lat'].values); clon = np.squeeze(d['lon'].values)
    C_up = np.squeeze(d['sfc_sw_up_all_clim'].values); C_dn = np.squeeze(d['sfc_sw_down_all_clim'].values)
with xr.open_dataset(SIC, decode_times=False) as d:
    sic = np.squeeze(d['sic'].values); slat = np.squeeze(d['latitude'].values)
if np.nanmax(sic) > 1.5: sic = sic / 100.0
if slat[0] > slat[-1]: sic = sic[:, ::-1, :]
sic = np.nan_to_num(sic, nan=0.0)                      # HadISST lon is already 0.5..359.5
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = lsm[0] if lsm.ndim == 3 else lsm
    mlat = np.squeeze(d['lat'].values); mlon = np.squeeze(d['lon'].values)
iy = np.abs(clat[:, None] - mlat[None, :]).argmin(1)
ix = np.abs(((clon[:, None] - mlon[None, :] + 180) % 360 - 180)).argmin(1)
rg = lambda a: a[..., iy, :][..., ix]
LAT = np.broadcast_to(clat[:, None], sic.shape[1:]); W = np.cos(np.deg2rad(LAT))
OCN = rg(lsm) <= 0.5

VARS = ('ci', 'ssr', 'ssrd', 'fal')


def load(arm):
    M = {v: 0.0 for v in VARS}
    for y in range(Y0, Y1 + 1):
        for v in VARS:
            with xr.open_dataset(f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc', decode_times=False) as d:
                k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
                M[v] = M[v] + rg(np.squeeze(d[k].values).astype('f8')) / (Y1 - Y0 + 1)
    return M


print(__doc__.split('Usage:')[0])
for arm in RUNS:
    M = load(arm)
    print(f'== {arm}  {Y0}-{Y1}   (both model and HadISST concentration >= {CONC})')
    print(f'  {"season":<22}{"Mkm2":>7}{"alb model":>11}{"fal":>7}{"alb CERES":>11}{"diff":>8}'
          f'{"SWdn mod":>9}{"SWdn CER":>9}{"absorbed":>10}{"abs CER":>9}{"per 0.01":>10}')
    for name, (mon, la, lb) in SEASON.items():
        ci = np.nan_to_num(np.stack([M['ci'][m] for m in mon]).mean(0))
        s = np.stack([sic[m] for m in mon]).mean(0)
        ssr = np.stack([M['ssr'][m] for m in mon]).mean(0)
        ssrd = np.stack([M['ssrd'][m] for m in mon]).mean(0)
        fal = np.stack([M['fal'][m] for m in mon]).mean(0)
        cup = np.stack([C_up[m] for m in mon]).mean(0); cdn = np.stack([C_dn[m] for m in mon]).mean(0)
        k = OCN & (ci >= CONC) & (s >= CONC) & (LAT >= la) & (LAT <= lb) & (ssrd > 1.0) & (cdn > 1.0)
        if k.sum() < 10:
            print(f'  {name:<22} fewer than 10 cells'); continue
        av = lambda f: float(np.average(f[k], weights=W[k]))
        am, ac = 1 - av(ssr) / av(ssrd), av(cup) / av(cdn)
        dn = av(ssrd) / 3600.0                            # accumulated J/m2 per output hour
        cd_ = av(cdn)                                      # CERES SWdn, already W/m2
        print(f'  {name:<22}{float((W * 111.195e3 ** 2)[k].sum() / 1e12):7.2f}{am:11.3f}{av(fal):7.3f}'
              f'{ac:11.3f}{am - ac:+8.3f}{dn:9.1f}{cd_:9.1f}{dn * (1 - am):10.1f}'
              f'{cd_ * (1 - ac):9.1f}{dn * 0.01:10.2f}')
print("""
  alb model = 1 - ssr/ssrd (flux-derived);  fal = the albedo OpenIFS used;  diff = model - CERES.
  SWdn and absorbed in W/m2 over the pack; "per 0.01" is what 0.01 of albedo is worth there.

  Namelist albedos in force (FESOM namelist.ice), against FESOM's own reference values
  commented in src/ice_modules.F90 and against in-situ literature:
    albsn  frozen snow   0.75   FESOM ref 0.81   observed dry snow on sea ice 0.80-0.87
    albsnm melting snow  0.65   FESOM ref 0.77   observed melting snow        0.70-0.77
    albi   frozen ice    0.66   FESOM ref 0.70   observed cold bare ice       0.60-0.70
    albim  melting ice   0.64   (no ref)         observed melting bare ice    0.50-0.60
    albpnd ponds         0.35                    observed ponds               0.15-0.35
  These four are POND-FREE surface albedos: ice_thermo_cpl.F90 builds alb_noponds from them and
  only then calls meltpond_albedo.  Darkening them to emulate ponds double-counts the ponds.
""")
