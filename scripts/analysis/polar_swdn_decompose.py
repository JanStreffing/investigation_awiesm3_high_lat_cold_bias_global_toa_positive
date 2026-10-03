"""Why is the polar pack short of downwelling shortwave: clear sky, or cloud?

ice_albedo_vs_ceres.py found the melt-season packs receive far less SWdn than CERES says
they should (NH JJA 161.6 against 212.4, SH DJF 212.7 against 253.0).  That deficit is six
times the energy a spring snowpack needs to melt, and it is the reason observationally
correct snow albedos send Arctic ice volume into a runaway: the dark albedos this campaign
used were compensating for it.

This splits the deficit.  SWdn_all = SWdn_clear - cloud attenuation, so comparing the model
and CERES clear-sky fields isolates whatever is wrong with the cloud-free atmosphere
(water vapour, aerosol, ozone, the solar geometry) from the cloud attenuation itself.  If
clear sky agrees and all sky does not, it is cloud, and the cloud columns then say whether
it is amount or opacity.

CAVEAT, the same one that applies to the albedo comparison: CERES surface fluxes over ice
are a radiative-transfer product rather than a measurement, and they are at their weakest
over bright surfaces.  CERES is also a 2005-2015 climatology against a pre-industrial run.
Read a 20 % discrepancy as real and a 5 % one as arguable.

Usage:  RUNS=PICAL_momixoff Y0=1940 Y1=1949 python3 scripts/analysis/polar_swdn_decompose.py
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('RUNS', 'PICAL_momixoff').split(',')[0]
Y0, Y1 = int(os.environ.get('Y0', 1940)), int(os.environ.get('Y1', 1949))
CONC = float(os.environ.get('CONC', 0.5))
CER = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
SIC = '/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc'
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'

with xr.open_dataset(CER, decode_times=False) as d:
    clat = np.squeeze(d['lat'].values); clon = np.squeeze(d['lon'].values)
    C = {k: np.squeeze(d[k].values) for k in
         ('sfc_sw_down_all_clim', 'sfc_sw_down_clr_t_clim', 'cldarea_total_daynight_clim',
          'cldtau_total_day_clim')}
with xr.open_dataset(SIC, decode_times=False) as d:
    sic = np.squeeze(d['sic'].values); slat = np.squeeze(d['latitude'].values)
if np.nanmax(sic) > 1.5: sic = sic / 100.0
if slat[0] > slat[-1]: sic = sic[:, ::-1, :]
sic = np.nan_to_num(sic, nan=0.0)
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = lsm[0] if lsm.ndim == 3 else lsm
    mlat = np.squeeze(d['lat'].values); mlon = np.squeeze(d['lon'].values)
iy = np.abs(clat[:, None] - mlat[None, :]).argmin(1)
ix = np.abs(((clon[:, None] - mlon[None, :] + 180) % 360 - 180)).argmin(1)
rg = lambda a: a[..., iy, :][..., ix]
LAT = np.broadcast_to(clat[:, None], sic.shape[1:]); W = np.cos(np.deg2rad(LAT))
OCN = rg(lsm) <= 0.5

VARS = ('ci', 'ssrd', 'ssrdc', 'tcc')
M = {}
for v in VARS:
    acc = 0.0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc', decode_times=False) as d:
            acc = acc + rg(np.squeeze(d[v].values).astype('f8')) / (Y1 - Y0 + 1)
    M[v] = acc

SEAS = {'NH melt (JJA)': ([5, 6, 7], 55, 90), 'SH melt (DJF)': ([11, 0, 1], -90, -55),
        'NH shoulder (MJ)': ([4, 5], 55, 90), 'SH shoulder (ON)': ([9, 10], -90, -55)}

print(__doc__.split('Usage:')[0])
print(f'{ARM}  {Y0}-{Y1}, pack cells where model and HadISST both >= {CONC}\n')
print(f'{"season":<19}{"SWdn clear":>21}{"SWdn all":>21}{"cloud attenuation":>23}{"cloud amount":>18}{"tau":>7}')
print(f'{"":<19}{"mod":>7}{"CER":>7}{"diff":>7}{"mod":>7}{"CER":>7}{"diff":>7}'
      f'{"mod":>8}{"CER":>8}{"diff":>7}{"mod":>9}{"CER":>9}{"CER":>7}')
for name, (mon, la, lb) in SEAS.items():
    ci = np.nan_to_num(np.stack([M['ci'][m] for m in mon]).mean(0))
    s = np.stack([sic[m] for m in mon]).mean(0)
    k = OCN & (ci >= CONC) & (s >= CONC) & (LAT >= la) & (LAT <= lb)
    if k.sum() < 10:
        print(f'  {name:<17} fewer than 10 cells'); continue
    av = lambda f: float(np.average(np.stack([f[m] for m in mon]).mean(0)[k], weights=W[k]))
    md_all, md_clr = av(M['ssrd']) / 3600.0, av(M['ssrdc']) / 3600.0
    cd_all, cd_clr = av(C['sfc_sw_down_all_clim']), av(C['sfc_sw_down_clr_t_clim'])
    ma, ca = md_clr - md_all, cd_clr - cd_all           # cloud attenuation of SWdn
    print(f'  {name:<17}{md_clr:7.1f}{cd_clr:7.1f}{md_clr-cd_clr:+7.1f}'
          f'{md_all:7.1f}{cd_all:7.1f}{md_all-cd_all:+7.1f}'
          f'{ma:8.1f}{ca:8.1f}{ma-ca:+7.1f}'
          f'{av(M["tcc"])*100:9.1f}{av(C["cldarea_total_daynight_clim"]):9.1f}'
          f'{av(C["cldtau_total_day_clim"]):7.1f}')
print("""
  Clear-sky diff near zero with a large all-sky diff means the cloud-free atmosphere is fine
  and the error is cloud.  The cloud-attenuation column is SWdn_clear - SWdn_all, i.e. how
  much shortwave the clouds remove; a positive diff there is the model removing too much.
  Cloud amount is in %; CERES tau is the total-cloud optical depth for reference.
""")
