"""Too much cloud water, or too many droplets?  Model against CLARA-A3 over the polar packs.

polar_swdn_decompose.py showed the model removes 50-57 W/m2 too much shortwave over the
Arctic melt-season pack while carrying roughly the observed cloud cover, and
polar_cloud_pack_vs_open.py showed the error belongs to the ice surface: at the same
latitude, open water is fine.  So the clouds over ice are too OPAQUE.  Optical thickness of
a liquid cloud goes as tau ~ LWP / r_eff, so the excess is either too much condensate or
droplets that are too small -- and those point at different levers.

CLARA-A3 (EUMETSAT CM SAF, AVHRR, /pool/data/ICDC) resolves that directly: it carries in-
cloud lwp, all-sky lwp, liquid optical thickness, effective radius AND droplet number, each
with an uncertainty estimate.

The model's droplet number is not prognostic: sucldp.F90 sets RCL_KK_CLOUD_NUM_SEA = 50 and
_LAND = 300 cm-3 as fixed constants, unset in fort.4, so the comparison against CLARA's
retrieved cdnc_liq is a direct test of a hard-coded number.

CAVEATS.  CLARA is a passive retrieval and cloud over snow and ice is its hardest case, the
same weakness that makes the CERES cloud fraction untrustworthy there; the uncertainty
columns are printed for that reason.  It needs daylight, so only the melt seasons are
usable.  And it is a present-day record against a pre-industrial run.

Usage:  RUNS=PICAL_momixoff Y0=1940 Y1=1949 python3 scripts/analysis/polar_lwp_vs_clara.py
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('RUNS', 'PICAL_momixoff').split(',')[0]
Y0, Y1 = int(os.environ.get('Y0', 1940)), int(os.environ.get('Y1', 1949))
CY0, CY1 = int(os.environ.get('CY0', 2008)), int(os.environ.get('CY1', 2014))
CONC = float(os.environ.get('CONC', 0.5))
CLARA = '/pool/data/ICDC/atmosphere/eumetsat_clara3_cloud/DATA/LiquidWaterPath/MONTHLY'
SIC = '/work/ab0246/a270092/obs/hadisst2/HadISST-SIC_monthly.nc'
LSMF = '/work/bb1469/a270092/runtime/awiesm3-v3.4/PI200/outdata/oifs/atm_remapped_1m_lsm_1390-1390.nc'
CVARS = ('lwp', 'lwp_allsky', 'lwp_unc_mean', 'cot_liq', 'cre_liq', 'cdnc_liq', 'nobs_liq_cot')


def clara_month(mm):
    """CY0..CY1 mean of month mm on the CLARA 0.5 deg grid."""
    acc, n = {v: 0.0 for v in CVARS}, 0
    for y in range(CY0, CY1 + 1):
        f = glob.glob(f'{CLARA}/LWPmm{y}{mm:02d}01*.nc')
        if not f:
            continue
        with xr.open_dataset(f[0], decode_times=False) as d:
            for v in CVARS:
                acc[v] = acc[v] + np.squeeze(d[v].values).astype('f8')
            global clat, clon
            clat, clon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
        n += 1
    if not n:
        raise SystemExit(f'no CLARA files for month {mm}')
    return {v: acc[v] / n for v in CVARS}


C = {mm: clara_month(mm) for mm in (6, 7, 8, 12, 1, 2)}

with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values); lsm = lsm[0] if lsm.ndim == 3 else lsm
    mlat = np.squeeze(d['lat'].values); mlon = np.squeeze(d['lon'].values)
with xr.open_dataset(SIC, decode_times=False) as d:
    sic = np.squeeze(d['sic'].values); slat = np.squeeze(d['latitude'].values)
    slon = np.squeeze(d['longitude'].values)
if np.nanmax(sic) > 1.5: sic = sic / 100.0
if slat[0] > slat[-1]: sic = sic[:, ::-1, :]
sic = np.nan_to_num(sic, nan=0.0)

# everything is evaluated on the CLARA grid; model and HadISST are expanded onto it by
# nearest index, the same convention the CERES scripts use
iy = np.abs(clat[:, None] - mlat[None, :]).argmin(1)
ix = np.abs(((clon[:, None] - mlon[None, :] + 180) % 360 - 180)).argmin(1)
rg = lambda a: a[..., iy, :][..., ix]
jy = np.abs(clat[:, None] - slat[None, :]).argmin(1)
jx = np.abs(((clon[:, None] - slon[None, :] + 180) % 360 - 180)).argmin(1)
rs = lambda a: a[..., jy, :][..., jx]

LAT = np.broadcast_to(clat[:, None], (clat.size, clon.size)); W = np.cos(np.deg2rad(LAT))
OCN = rg(lsm) <= 0.5

M = {}
for v in ('ci', 'tclw', 'tciw', 'lcc', 'tcc'):
    a = 0.0
    for y in range(Y0, Y1 + 1):
        with xr.open_dataset(f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc', decode_times=False) as d:
            a = a + rg(np.squeeze(d[v].values).astype('f8')) / (Y1 - Y0 + 1)
    M[v] = a

print(__doc__.split('Usage:')[0])
print(f'model {ARM} {Y0}-{Y1}   CLARA-A3 {CY0}-{CY1}   pack: model and HadISST both >= {CONC}\n')
hdr = (f'{"season":<17}{"LWP grid-mean":>20}{"LWP in-cloud":>21}{"r_eff":>15}{"CDNC":>16}')
print(hdr)
print(f'{"":<17}{"mod":>8}{"CLARA":>7}{"unc":>5}{"mod":>9}{"CLARA":>7}{"tau":>5}'
      f'{"CLARA":>8}{"":>7}{"mod":>8}{"CLARA":>8}')
print(f'{"":<17}{"g/m2":>8}{"g/m2":>7}{"g/m2":>5}{"g/m2":>9}{"g/m2":>7}{"":>5}'
      f'{"um":>8}{"":>7}{"cm-3":>8}{"cm-3":>8}')
for name, (mon, la, lb) in {'NH melt (JJA)': ([6, 7, 8], 55, 90),
                            'SH melt (DJF)': ([12, 1, 2], -90, -55)}.items():
    mi = [m - 1 for m in mon]                                  # model arrays are 0-based
    ci_ = np.nan_to_num(np.stack([M['ci'][m] for m in mi]).mean(0))
    s = np.stack([rs(sic[m]) for m in mi]).mean(0)
    nob = np.stack([C[m]['nobs_liq_cot'] for m in mon]).mean(0)
    k = OCN & (ci_ >= CONC) & (s >= CONC) & (LAT >= la) & (LAT <= lb) & np.isfinite(nob) & (nob >= 3)
    if k.sum() < 10:
        print(f'  {name:<15} fewer than 10 cells with CLARA sampling'); continue
    av = lambda f: float(np.average(np.nan_to_num(f)[k], weights=W[k]))
    cav = lambda v: av(np.stack([C[m][v] for m in mon]).mean(0))
    m_lwp = av(np.stack([M['tclw'][m] for m in mi]).mean(0)) * 1e3          # kg/m2 -> g/m2
    m_lcc = av(np.stack([M['lcc'][m] for m in mi]).mean(0))
    print(f'  {name:<15}{m_lwp:8.1f}{cav("lwp_allsky")*1e3:7.1f}{cav("lwp_unc_mean")*1e3:5.1f}'
          f'{m_lwp/max(m_lcc,1e-3):9.1f}{cav("lwp")*1e3:7.1f}{cav("cot_liq"):5.1f}'
          f'{cav("cre_liq")*1e6:8.1f}{"":>7}{50.0:8.1f}{cav("cdnc_liq")/1e6:8.1f}')
print("""
  Model LWP is tclw, the grid-mean total column liquid; in-cloud divides it by low cloud
  cover.  Model CDNC is not retrieved, it is the hard-coded RCL_KK_CLOUD_NUM_SEA = 50 cm-3
  (sucldp.F90:456), shown for comparison with what CLARA sees.
  Too much LWP points at condensate or precipitation efficiency; too small r_eff / too high
  CDNC points at the droplet number, which over sea ice is a constant nobody has revisited.
""")
