"""Is the polar cloud excess tied to the ICE SURFACE, or is it the cloud scheme?

polar_swdn_decompose.py found the model puts 94.3 % cloud over the Arctic melt-season pack
against CERES's 80.6 %, while the report's standing table has the model 6.7 points SHORT of
cloud at 45-65S.  A global cloud-fraction knob therefore fixes one and worsens the other.

A selective fix exists only if the Arctic excess belongs to the ice surface rather than to
high latitudes in general.  This compares cloud cover over pack and over open water WITHIN
THE SAME LATITUDE BAND, so insolation, solar zenith angle and the large-scale flow are held
as fixed as they can be, and the only difference left is what is underneath the cloud.

  pack  = model ci >= 0.5  and HadISST sic >= 0.5   (both agree there is ice)
  open  = model ci <= 0.15 and HadISST sic <= 0.15  (both agree there is not)

If the model-minus-CERES gap is large over pack and small over open water, the error is in
the ice-atmosphere coupling and can be attacked without touching the Southern Ocean.  If it
is the same over both, it is the cloud scheme and any fix is a trade.

Usage:  RUNS=PICAL_momixoff Y0=1940 Y1=1949 python3 scripts/analysis/polar_cloud_pack_vs_open.py
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('RUNS', 'PICAL_momixoff').split(',')[0]
Y0, Y1 = int(os.environ.get('Y0', 1940)), int(os.environ.get('Y1', 1949))
CER = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
SIC = '/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc'
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'

with xr.open_dataset(CER, decode_times=False) as d:
    clat = np.squeeze(d['lat'].values); clon = np.squeeze(d['lon'].values)
    C = {k: np.squeeze(d[k].values) for k in
         ('cldarea_total_daynight_clim', 'sfc_sw_down_all_clim', 'sfc_sw_down_clr_t_clim')}
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

M = {}
for v in ('ci', 'tcc', 'ssrd', 'ssrdc'):
    a = 0.0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc', decode_times=False) as d:
            a = a + rg(np.squeeze(d[v].values).astype('f8')) / (Y1 - Y0 + 1)
    M[v] = a

BANDS = [('NH JJA', [5, 6, 7], [(55, 65), (65, 75), (75, 90)]),
         ('SH DJF', [11, 0, 1], [(-65, -55), (-75, -65), (-90, -75)])]

print(__doc__.split('Usage:')[0])
print(f'{ARM}  {Y0}-{Y1}\n')
print(f'{"band":<16}{"surface":<7}{"cells":>7}{"cloud mod":>10}{"cloud CER":>10}{"diff":>7}'
      f'{"atten mod":>10}{"atten CER":>10}{"diff":>7}')
for season, mon, bands in BANDS:
    ci = np.nan_to_num(np.stack([M['ci'][m] for m in mon]).mean(0))
    s = np.stack([sic[m] for m in mon]).mean(0)
    for la, lb in bands:
        inband = OCN & (LAT >= la) & (LAT <= lb)
        for lab, k in (('pack', inband & (ci >= 0.5) & (s >= 0.5)),
                       ('open', inband & (ci <= 0.15) & (s <= 0.15))):
            if k.sum() < 10:
                print(f'  {season} {la:>3}/{lb:<4}{lab:<7}{int(k.sum()):>7}   fewer than 10 cells'); continue
            av = lambda f: float(np.average(np.stack([f[m] for m in mon]).mean(0)[k], weights=W[k]))
            mc, cc = av(M['tcc']) * 100, av(C['cldarea_total_daynight_clim'])
            ma = (av(M['ssrdc']) - av(M['ssrd'])) / 3600.0
            ca = av(C['sfc_sw_down_clr_t_clim']) - av(C['sfc_sw_down_all_clim'])
            print(f'  {season} {la:>3}/{lb:<4}{lab:<7}{int(k.sum()):>7}{mc:10.1f}{cc:10.1f}{mc-cc:+7.1f}'
                  f'{ma:10.1f}{ca:10.1f}{ma-ca:+7.1f}')
print("""
  cloud is % cover; atten is SWdn_clear - SWdn_all in W/m2, i.e. what the cloud removes.
  A gap that is large over pack and small over open water at the SAME latitude points at the
  ice surface, not the cloud scheme.
""")
